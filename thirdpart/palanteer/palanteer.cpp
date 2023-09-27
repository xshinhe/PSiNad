#define PL_IMPLEMENTATION                  1
#include "palanteer.h"

//=========================================================================

#if USE_PL==1

namespace plPriv {

GlobalContext_t globalCtx(PL_IMPL_DYN_STRING_QTY);

thread_local ThreadContext_t threadCtx;

#if PL_NOASSERT==0

void
registerCli(plCliHandler_t handler, const char* name, const char* specParams, const char* description,
            hashStr_t nameHash, hashStr_t specParamsHash, hashStr_t descriptionHash)
{
    implCtx.cliManager.registerCli(handler, name, specParams, description,
                                   nameHash, specParamsHash, descriptionHash);
}

static uint8_t*
helperFillResponseBufferHeader(RemoteCommandType commandType, int commandByteSize, uint8_t* br)
{
    br[0] = 'P';
    br[1] = 'L';
    br[2] = ((int)PL_DATA_TYPE_CONTROL>>8)&0xFF;
    br[3] = ((int)PL_DATA_TYPE_CONTROL>>0)&0xFF;
    commandByteSize += 2;               // Size of the command type
    br[4] = (commandByteSize>>24)&0xFF; // Command byte quantity (after the 8 bytes remote data type header)
    br[5] = (commandByteSize>>16)&0xFF;
    br[6] = (commandByteSize>> 8)&0xFF;
    br[7] = (commandByteSize>> 0)&0xFF;
    br[8] = ((int)commandType>>8)&0xFF;
    br[9] = ((int)commandType>>0)&0xFF;
    return br;
}

static void
helperFinishResponseBuffer(int txBufferSize)
{
    // Mark it for sending (in the tx thread)
    std::lock_guard<std::mutex> lk(implCtx.txThreadSyncMx);
    implCtx.rspBufferSize.store(txBufferSize);
    // Will be sent at the next tx thread loop, that we force now
    implCtx.txThreadSyncCv.notify_one();
}

void
receiveFromServer(void)
{
    auto& ic = implCtx;
    plgDeclareThread(PL_VERBOSE, "Palanteer/Reception");

    while(!ic.threadServerFlagStop.load()) {

        int recByteQty = palComReceive(ic.reqBuffer, PL_IMPL_REMOTE_REQUEST_BUFFER_BYTE_QTY);
        if     (recByteQty <0) continue; // Timeout on reception (empty)
        else if(recByteQty==0) break;    // Client is disconnected

        // Parse the received content, expected block structure is:
        //   [block]            <2B: synchro magic>: 'P' 'L'
        //   [block]            <2B: bloc type>
        //   [remote data type] <4B: command byte qty>
        //   [remote data type] <2B: remote command type>
        //   (following payloads depends on the command type)
        // Sanity
        uint8_t* b = (uint8_t*)ic.reqBuffer;
        plAssert(recByteQty>=10);
        plAssert(b[0]=='P' && b[1]=='L', "Magic not present: connection is desynchronized");
        DataType dt = (DataType)((b[2]<<8) | b[3]);
        plAssert(dt==PL_DATA_TYPE_CONTROL, "Wrong block data type received through socket despite synchronization: connection is buggy.");
        int commandByteQty = (b[4]<<24) | (b[5]<<16) | (b[6]<<8) | b[7];
        plAssert(8+commandByteQty<=PL_IMPL_REMOTE_REQUEST_BUFFER_BYTE_QTY, "Too big remote command received. Limit is:",
                 PL_IMPL_REMOTE_REQUEST_BUFFER_BYTE_QTY, 8+commandByteQty);

        // Wait the end of the request reception
        while(recByteQty<8+commandByteQty && !ic.threadServerFlagStop.load()) {
            int nextByteQty = palComReceive(ic.reqBuffer+recByteQty, PL_IMPL_REMOTE_REQUEST_BUFFER_BYTE_QTY-recByteQty);
            if     (nextByteQty <0) continue; // Timeout on reception (empty)
            else if(nextByteQty==0) break;    // Client is disconnected
            recByteQty += nextByteQty;
        }
        if(recByteQty<8+commandByteQty || ic.threadServerFlagStop.load()) continue; // Exit or connection break case
        if(ic.rspBufferSize.load()>0) continue; // The response buffer shall be free, if the sender behaves as expected

        // Process the command
        RemoteCommandType ct = (RemoteCommandType)((b[8]<<8) | b[9]);
        int payloadByteQty = commandByteQty-2; // 2 = command type

        if(ct==PL_CMD_SET_FREEZE_MODE) {
            { // Scope so that the bootstrap of the thread is not included
                plgScope(PL_VERBOSE, "Request: set freeze mode");
                plAssert(payloadByteQty==1);

                // Update the state
                bool isFreezeEnabled = (b[10]!=0);
                plgData(PL_VERBOSE, "State", isFreezeEnabled);
                {
                    std::unique_lock<std::mutex> lk(ic.frozenThreadMx);
                    ic.frozenThreadEnabled.store(isFreezeEnabled? 1:0);
                }

                // Free the threads if disabled
                if(!isFreezeEnabled) {
                    uint64_t bitmap = ic.frozenThreadBitmap.load();
                    int tId = 0;
                    while(bitmap) {
                        if(bitmap&1) ic.frozenThreadCv[tId].notify_one();
                        bitmap >>= 1; ++tId;
                    }
                }

                // Build and send the response
                uint8_t* br = helperFillResponseBufferHeader(PL_CMD_SET_FREEZE_MODE, 2, ic.rspBuffer);
                br[10] = (((int)PL_OK)>>8)&0xFF;
                br[11] = (((int)PL_OK)>>0)&0xFF;
                helperFinishResponseBuffer(12);
            }

            // Notify the initialization thread that server and reception thread is ready
            // This let also the chance to safely activate the freeze mode before starting the program
            if(!ic.rxIsStarted) {
                std::lock_guard<std::mutex> lk(ic.threadInitMx);
                ic.rxIsStarted = true;
                ic.threadInitCvTx.notify_one();
            }
        }

        else if(ct==PL_CMD_STEP_CONTINUE) {
            plgScope(PL_VERBOSE, "Request: resume thread execution");
            plAssert(payloadByteQty==8);
            // Unmask the selected threads
            uint64_t bitmap = ((uint64_t)b[10]<<56) | ((uint64_t)b[11]<<48) | ((uint64_t)b[12]<<40) | ((uint64_t)b[13]<<32) |
                ((uint64_t)b[14]<<24) | ((uint64_t)b[15]<<16) | ((uint64_t)b[16]<<8) | ((uint64_t)b[17]<<0);
            plgData(PL_VERBOSE, "Thread bitmap##hexa", bitmap);
            {
                std::unique_lock<std::mutex> lk(ic.frozenThreadMx);
                ic.frozenThreadBitmap.fetch_and(~bitmap);
            }
            // Wake up the selected threads
            int tId = 0;
            while(bitmap) {
                if(bitmap&1) ic.frozenThreadCv[tId].notify_one();
                bitmap >>= 1; ++tId;
            }
            // Build and send the response
            uint8_t* br = helperFillResponseBufferHeader(PL_CMD_STEP_CONTINUE, 2, ic.rspBuffer);
            br[10] = (((int)PL_OK)>>8)&0xFF;
            br[11] = (((int)PL_OK)>>0)&0xFF;
            helperFinishResponseBuffer(12);
        }

        else if(ct==PL_CMD_SET_MAX_LATENCY) {
            plgScope(PL_VERBOSE, "Request: set max latency");
            plAssert(payloadByteQty==2);
            // Update the internal state
            int maxLatencyMs = (int)(b[10]<<8) | (int)b[11];
            plgData(PL_VERBOSE, "Max latency##ms", maxLatencyMs);
            ic.maxSendingLatencyNs = 1e6*maxLatencyMs;

            // Build and send the response
            uint8_t* br = helperFillResponseBufferHeader(PL_CMD_SET_MAX_LATENCY, 2, ic.rspBuffer);
            br[10] = (((int)PL_OK)>>8)&0xFF;
            br[11] = (((int)PL_OK)>>0)&0xFF;
            helperFinishResponseBuffer(12);
        }

        else if(ct==PL_CMD_KILL_PROGRAM) {
            plgScope(PL_VERBOSE, "Request: kill program");
            // We use quick_exit (introduced in C++11) instead of exit (or abort) because:
            //  - 'quick_exit' is design for multi-threaded application which are hard or costly in plumbing to stop in a clean manner
            //  - 'quick_exit' does not "clean" the process before leaving but allows manual cleaning through 'at_quick_exit' registration
            //  - 'exit' has good odds to block or crash if the program is not specifically designed for it
            //  - 'abort' is "violent", the signal does not allow any cleaning and may lead to a core dump or a popup window (on windows)
            std::quick_exit(0);
            // No need to bother for any response, as connection will be down very soon...
        }

        else if(ct==PL_CMD_CALL_CLI) {
            plgScope(PL_VERBOSE, "Request: CLI");
            // Sanity
            plAssert(payloadByteQty>2);
            b[8+commandByteQty-1] = 0; // Force the zero terminated string at the end of the reception buffer, just in case

            // Get the CLI request qty and prepare the response buffer
            int cliRequestQty = (b[10]<<8) | b[11];
            int reqOffset     = 12; // Current position in request buffer
            int rspOffset     = 12; // Command response header, partially filled at the end
            uint8_t* br       = helperFillResponseBufferHeader(PL_CMD_CALL_CLI, 0, ic.rspBuffer); // 0 bytes size is a dummy value, it will be overwritten later
            // br[ 4-> 7] = 4 bytes for the command response byte qty (filled later)
            // br[10->11] = 2 bytes for the CLI response qty (filled later)
            plgData(PL_VERBOSE, "CLI request quantity", cliRequestQty);

            // Loop on requests and fill the response buffer
            constexpr int bufferFullMessageLength = 28; //strlen("CLI response buffer is full")+1; // strlen is not constexpr in MSVC2019
            int cliRequestNbr = 0;
            while(cliRequestNbr<cliRequestQty) {
                // Call
                plRemoteStatus cliStatus;
                plgScope(PL_VERBOSE, "Call");
                const char* cliResponse = ic.cliManager.execute((char*)&b[reqOffset], cliStatus);
                int responseLength = (int)strlen(cliResponse)+1;
                plgVar(PL_VERBOSE, cliRequestNbr, cliStatus, responseLength);

                // Check if we can store the response in the buffer
                if(rspOffset+2+responseLength>PL_IMPL_REMOTE_RESPONSE_BUFFER_BYTE_QTY -
                   ((cliRequestNbr==cliRequestQty-1)? 0 : 2+bufferFullMessageLength)) { // minimum size to store a truncated response, if not last command
                    // Not enough space in response buffer
                    plAssert(rspOffset+2+bufferFullMessageLength<=PL_IMPL_REMOTE_RESPONSE_BUFFER_BYTE_QTY); // It should by design
                    plgLogError(PL_VERBOSE, "error", "Not enough space in the response buffer");
                    br[rspOffset+0] = (((int)PL_ERROR)>>8)&0xFF;
                    br[rspOffset+1] = (((int)PL_ERROR)>>0)&0xFF;
                    snprintf((char*)&br[rspOffset+2], bufferFullMessageLength+1, "CLI response buffer is full");
                    rspOffset += 2+bufferFullMessageLength;
                    while(b[reqOffset]!=0) ++reqOffset;
                    ++reqOffset; // Skip the zero termination of the string request
                    ++cliRequestNbr;
                    break;
                }

                // Store the answer
                br[rspOffset+0] = (((int)cliStatus)>>8)&0xFF;
                br[rspOffset+1] = (((int)cliStatus)>>0)&0xFF;
                memcpy(br+rspOffset+2, cliResponse, responseLength); // Copy the response string with the zero termination
                rspOffset += 2+responseLength;

                // Next request
                while(b[reqOffset]!=0) ++reqOffset;
                ++reqOffset; // Skip the zero termination of the string request
                ++cliRequestNbr;

                // if(cliStatus!=PL_OK) break;  // Another possible behavior is stop at first failure (not by default)
            } // End of loop on the CLI requests
            // Either not all requests are processed, either we read all the request bytes
            plAssert(cliRequestNbr<cliRequestQty || reqOffset==8+commandByteQty,
                     cliRequestNbr, cliRequestQty, reqOffset, 8+commandByteQty);

            // Finalize the response and send it
            br[4] = ((rspOffset-8)>>24)&0xFF; // Command byte quantity (after the 8 bytes remote data type header)
            br[5] = ((rspOffset-8)>>16)&0xFF;
            br[6] = ((rspOffset-8)>> 8)&0xFF;
            br[7] = ((rspOffset-8)>> 0)&0xFF;
            br[10] = (cliRequestNbr>>8)&0xFF; // CLI answer quantity (less than or equal to the request quantity)
            br[11] = (cliRequestNbr>>0)&0xFF;
            helperFinishResponseBuffer(rspOffset);
        } // if(ct==PL_CMD_CALL_CLI)

    } // End of reception loop

    // In case of server connection failure, the program shall be started anyway
    if(!ic.rxIsStarted) {
        std::lock_guard<std::mutex> lk(ic.threadInitMx);
        ic.rxIsStarted = true;
        ic.threadInitCvTx.notify_one();
    }

    plgLogInfo(PL_VERBOSE, "threading", "End of Palanteer reception loop");
}



#if PL_EXTERNAL_STRINGS==0
    void PL_NOINLINE
    failedAssertSimple(const char* filename, int lineNbr, const char* function, const char* condition)
    {
        char infoStr[1024];
        snprintf(infoStr, sizeof(infoStr), "[PALANTEER] Assertion failed: %s\n  On function: %s\n  On file    : %s(%d)\n", condition, function, filename, lineNbr);
        plCrash(infoStr);
    }
#else
    void PL_NOINLINE
    failedAssertSimpleEs(hashStr_t filenameHash, int lineNbr, hashStr_t conditionHash)
    {
        char infoStr[1024];
        snprintf(infoStr, sizeof(infoStr), "[PALANTEER] Assertion failed: @@%016" PL_PRI_HASH "@@\n  On file @@%016" PL_PRI_HASH "@@(%d)\n",
                 conditionHash, filenameHash, lineNbr);
        plCrash(infoStr);
    }
#endif // if PL_EXTERNAL_STRINGS==0
#endif // if PL_NOASSERT==0


} // namespace plPriv


// @#LATER Handle the alignments stuff
void* operator new  (std::size_t size) noexcept(false)                 { void* ptr = PL_NEW_(ptr, size); return(ptr); }
void* operator new[](std::size_t size) noexcept(false)                 { void* ptr = PL_NEW_(ptr, size); return(ptr); }
void* operator new  (std::size_t size, const std::nothrow_t &) throw() { void* ptr = PL_NEW_(ptr, size); return(ptr); }
void* operator new[](std::size_t size, const std::nothrow_t &) throw() { void* ptr = PL_NEW_(ptr, size); return(ptr); }

void operator delete  (void* ptr) noexcept                        { PL_DELETE_(ptr); }
void operator delete[](void* ptr) noexcept                        { PL_DELETE_(ptr); }
void operator delete  (void* ptr, std::size_t size) noexcept      { PL_DELETE_(ptr); PL_UNUSED(size); }
void operator delete[](void* ptr, std::size_t size) noexcept      { PL_DELETE_(ptr); PL_UNUSED(size); }
void operator delete  (void *ptr, const std::nothrow_t&) noexcept { PL_DELETE_(ptr); }
void operator delete[](void *ptr, const std::nothrow_t&) noexcept { PL_DELETE_(ptr); }

#if PL_NOCONTROL==0
void
plFreezePoint(void)
{
    // Is freeze mode enabled?
    auto& ic = plPriv::implCtx;
    if(!ic.frozenThreadEnabled.load()) return;

    // Mark the thread as frozen
    int tId = plPriv::getThreadId();
    if(tId>=PL_MAX_THREAD_QTY) return; // No freeze feature above the maximum thread quantity
    uint64_t mask = 1ULL<<tId;
    ic.frozenThreadBitmap.fetch_or(mask);
    ic.frozenThreadBitmapChange.fetch_or(mask);

    // Wait for unfreeze
    std::unique_lock<std::mutex> lk(ic.frozenThreadMx);
    ic.frozenThreadCv[tId].wait(lk, [&] { return !(ic.frozenThreadEnabled.load() && (ic.frozenThreadBitmap.load()&mask)); });
    // Unfrozen, force the thread bit to zero in case the freeze feature has been disabled
    ic.frozenThreadBitmap.fetch_and(~mask);
    ic.frozenThreadBitmapChange.fetch_or(mask);
}
#endif // if PL_NOCONTROL==0

void
plSetFilename(const char* filename)
{
    snprintf(plPriv::implCtx.filename, sizeof(plPriv::implCtx.filename), "%s", filename); // Safe and null terminated
}


void
plSetServer(const char* serverAddr, int serverPort)
{
    snprintf(plPriv::implCtx.serverAddr, sizeof(plPriv::implCtx.serverAddr), "%s", serverAddr);
    plPriv::implCtx.serverPort = serverPort;
}


void
plSetLogLevelRecord(plLogLevel level)
{
    plAssert(level>=PL_LOG_LEVEL_ALL && level<=PL_LOG_LEVEL_NONE, level);
    plPriv::globalCtx.minLogLevelRecord = level;
}


void
plSetLogLevelConsole(plLogLevel level)
{
    plAssert(level>=PL_LOG_LEVEL_ALL && level<=PL_LOG_LEVEL_NONE, level);
    plPriv::globalCtx.minLogLevelConsole = level;
}


void
plInitAndStart(const char* appName, plMode mode, const char* buildName, int serverConnectionTimeoutMsec)
{
    auto& ic = plPriv::implCtx;
    ic.mode = mode;
    (void)buildName;

    // Sanity
    static_assert(PL_MAX_THREAD_QTY<=254, "Maximum supported thread quantity reached (limitation on exchange structure side)");
    static_assert(PL_IMPL_COLLECTION_BUFFER_BYTE_QTY>(int)2*sizeof(plPriv::EventInt), "Too small collection buffer"); // Much more expected anyway...
    static_assert(PL_IMPL_DYN_STRING_QTY>=32, "Invalid configuration");  // Stack trace requires dynamic strings
#if PL_NOCONTROL==0 || PL_NOEVENT==0
#if PL_COMPACT_MODEL==1
    static_assert(sizeof(plPriv::EventExt)==12, "Bad size of compact exchange event structure");
#else
    static_assert(sizeof(plPriv::EventExt)==24, "Bad size of full exchange event structure");
#endif
    plAssert(!ic.allocCollectBuffer, "Double call to 'plInitAndStart' detected");
#if PL_IMPL_AUTO_INSTRUMENT==1 && (!defined(__GNUC__) || defined(__clang__))
    // Clang misses the file exclusion feature (-finstrument-functions-exclude-file-list) which is mandatory to avoid
    // "Larsen effect" by excluding palanteer.h and all functions from the standard library used when logging (atomic, thread...).
    #error "Sorry, auto instrumentation is supported only for GCC"
#endif
#endif

    // Register POSIX signals
    memset(ic.signalsOldHandlers, 0, sizeof(ic.signalsOldHandlers));
#if PL_IMPL_CATCH_SIGNALS==1
    ic.signalsOldHandlers[0] = std::signal(SIGABRT, plPriv::signalHandler);
    ic.signalsOldHandlers[1] = std::signal(SIGFPE,  plPriv::signalHandler);
    ic.signalsOldHandlers[2] = std::signal(SIGILL,  plPriv::signalHandler);
    ic.signalsOldHandlers[3] = std::signal(SIGSEGV, plPriv::signalHandler);
    ic.signalsOldHandlers[4] = std::signal(SIGINT,  plPriv::signalHandler);
    ic.signalsOldHandlers[5] = std::signal(SIGTERM, plPriv::signalHandler);
#if defined(__unix__)
    ic.signalsOldHandlers[6] = std::signal(SIGPIPE, plPriv::signalHandler);
#endif
    ic.signalHandlersSaved = true;
#if defined(_WIN32)
    // Register the exception handler
    ic.exceptionHandler = AddVectoredExceptionHandler(1, plPriv::exceptionHandler);
#endif // if defined(_WIN32)
#endif // if PL_IMPL_CATCH_SIGNALS==1

#if defined(_WIN32) && PL_IMPL_STACKTRACE==1
    // Initialize the symbol reading for the stacktrace (in case of crash)
    SymInitialize(GetCurrentProcess(), 0, true);
    SymSetOptions(SYMOPT_LOAD_LINES);
    ic.rtlWalkFrameChain = (plPriv::rtlWalkFrameChain_t)GetProcAddress(GetModuleHandleA("ntdll.dll"), "RtlWalkFrameChain");
    plAssert(ic.rtlWalkFrameChain);
#endif // if defined(_WIN32) && PL_IMPL_STACKTRACE==1

    // In inactive mode, nothing else to do
    if(mode==PL_MODE_INACTIVE) {
        return;
    }

    // Remove some warning for some configurations
    PL_UNUSED(appName);
    PL_UNUSED(serverConnectionTimeoutMsec);

#if PL_NOEVENT==0 && PL_VIRTUAL_THREADS==1
    for(int i=0; i<PL_MAX_THREAD_QTY; ++i) {
        // Other fields are not reset to have persistent thread names
        plPriv::globalCtx.threadInfos[i].isSuspended = false;
        plPriv::globalCtx.threadInfos[i].isBeginSent = false;
    }
#endif

#if PL_IMPL_AUTO_INSTRUMENT==1 && defined(__unix__) && (PL_NOCONTROL==0 || PL_NOEVENT==0)
    // Initialize the automatic instrumentation mode (gcc only) with -finstrument-functions
    ic.hashedAppName = plPriv::hashString(appName);
#if PL_IMPL_AUTO_INSTRUMENT_PIC==1
    Dl_info info;
    if(dladdr((void*)&plInitAndStart, &info)) {
        ic.baseVirtualAddress = (char*)info.dli_fbase;
    }
#else
    ic.baseVirtualAddress = 0;  // Absolute addresses
#endif  // PL_IMPL_AUTO_INSTRUMENT_PIC==1
#endif  // PL_IMPL_AUTO_INSTRUMENT==1 && defined(__unix__) && (PL_NOCONTROL==0 || PL_NOEVENT==0)

#if PL_NOCONTROL==0 || PL_NOEVENT==0
    // Measure the event's high performance clock frequency with the standard nanosecond clock
    int64_t highPerfT0 = PL_GET_CLOCK_TICK_FUNC();
    const auto stdT0 = std::chrono::high_resolution_clock::now();
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    int64_t highPerfT1 = PL_GET_CLOCK_TICK_FUNC();
    const auto stdT1 = std::chrono::high_resolution_clock::now();
    // This coefficient will be sent to the server
#if PL_SHORT_DATE==1
    if(highPerfT1<=highPerfT0) highPerfT1 += (1LL<<32);
#endif
    ic.tickToNs = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(stdT1-stdT0).count()/(double)(highPerfT1-highPerfT0);

#if defined(_WIN32) && PL_NOEVENT==0 && PL_IMPL_CONTEXT_SWITCH==1
    // Windows: compute the clock conversion information for the context switch
    LARGE_INTEGER qpc;
    QueryPerformanceCounter(&qpc);
    ic.rdtscRef = plPriv::getClockTick();
    ic.qpcRef   = qpc.QuadPart;
    QueryPerformanceFrequency(&qpc);
    ic.qpcToRdtsc = 1e9/(qpc.QuadPart*ic.tickToNs);
#endif

    // Allocate the 2 collection banks (in one chunk, with a slight shift for a more efficient collectEvents())
    //   aligned on 64 bytes to match most cache lines (the internal event representation has also a size of 64 bytes)
    plPriv::globalCtx.collectBufferMaxEventQty = PL_IMPL_COLLECTION_BUFFER_BYTE_QTY/sizeof(plPriv::EventInt);
#if PL_NOEVENT==0
    plAssert((uint32_t)plPriv::globalCtx.collectBufferMaxEventQty<plPriv::EVTBUFFER_MASK_INDEX, "The collection buffer is too large");
#endif
    const int realBufferEventQty = plPriv::globalCtx.collectBufferMaxEventQty + (1+PL_MAX_THREAD_QTY)+64; // 64=margin for the collection thread
    int sendBufferSize = sizeof(plPriv::EventExt)*realBufferEventQty;
    if(sendBufferSize<PL_IMPL_REMOTE_RESPONSE_BUFFER_BYTE_QTY) sendBufferSize = PL_IMPL_REMOTE_RESPONSE_BUFFER_BYTE_QTY;
    ic.sendBuffer = new uint8_t[sendBufferSize+64];  // 64 = sent header margin
    ic.allocCollectBuffer = new uint8_t[sizeof(plPriv::EventInt)*2*realBufferEventQty+64];
    memset(ic.allocCollectBuffer, 0, sizeof(plPriv::EventInt)*2*realBufferEventQty+64);
    uint8_t* alignedAllocCollectBuffer = (uint8_t*)((((uintptr_t)ic.allocCollectBuffer)+64)&(uintptr_t)(~0x3F));
    plPriv::globalCtx.collectBuffers[0] = (plPriv::EventInt*)alignedAllocCollectBuffer;
    plPriv::globalCtx.collectBuffers[1] = plPriv::globalCtx.collectBuffers[0] + realBufferEventQty;

    // Initialize some fields
    memset(&ic.stats, 0, sizeof(plStats));
    ic.stats.collectBufferSizeByteQty = PL_IMPL_COLLECTION_BUFFER_BYTE_QTY;
    ic.stats.collectDynStringQty      = PL_IMPL_DYN_STRING_QTY;
    plPriv::globalCtx.bankAndIndex.store(0);
    ic.prevBankAndIndex = 1UL<<31;

    plPriv::palComInit(serverConnectionTimeoutMsec);
    if(ic.mode==PL_MODE_INACTIVE) return;


#if defined(__unix__) && PL_NOEVENT==0 && PL_IMPL_CONTEXT_SWITCH==1
    {
        ic.cswitchPollEnabled = true;
        // Configure the tracing
        char tmpStr[64];
        int  tracerFd;
        PL_WRITE_TRACE_("tracing_on",     "0",                   true); // Disables tracing while configuring
        PL_WRITE_TRACE_("current_tracer", "nop",                 true); // Removes all function tracers
        PL_WRITE_TRACE_("trace_options",  "noirq-info",         false); // No need for irq information
        PL_WRITE_TRACE_("trace_options",  "noannotate",         false); // No need for extra "annotate" infos
        PL_WRITE_TRACE_("trace_options",  "norecord-cmd",       false); // No need for extra thread info (PID are enough for our usage)
        PL_WRITE_TRACE_("trace_options",  "norecord-tgid",      false); // No need for the Thread Group ID
#if defined(__x86_64__)
        PL_WRITE_TRACE_("trace_clock",    "x86-tsc",             true); // Same clock than the default RDTSC one for Linux
#else
        PL_WRITE_TRACE_("trace_clock",    "mono",                true); // Usually (depending on arch) same as std::chrono::steady_clock but lower precision than RDTSC
#endif
        PL_WRITE_TRACE_("events/enable",  "0",                  false); // Disable all kernel events
        PL_WRITE_TRACE_("events/sched/sched_switch/enable", "1", true); // Enable the events we want
        PL_WRITE_TRACE_("events/irq/softirq_entry/enable",  "1", true);
        PL_WRITE_TRACE_("events/irq/softirq_exit/enable",   "1", true);
        PL_WRITE_TRACE_("buffer_size_kb", "512",                 true); // Reserve 512KB for exchanges
        PL_WRITE_TRACE_("tracing_on",     "1",                   true); // Enable tracing

        // Open the exchange pipe
        if(ic.cswitchPollEnabled && (ic.cswitchPollFd.fd=open("/sys/kernel/debug/tracing/trace_pipe", O_RDONLY))>=0) {
            ic.cswitchPollFd.events = POLLIN | POLLERR;
            ic.cswitchPollBuffer    = new char[plPriv::SWITCH_CTX_BUFFER_SIZE];
        } else ic.cswitchPollEnabled = false;
    }
#endif //if defined(__unix__) && PL_NOEVENT==0 && PL_IMPL_CONTEXT_SWITCH==1

#if defined(_WIN32) && PL_NOEVENT==0 && PL_IMPL_CONTEXT_SWITCH==1
    // See https://caseymuratori.com/blog_0025: "The Worst API Ever Made"

    // Allocate the tracer "property" structure as intended by the API
    ULONG propertySize = sizeof(EVENT_TRACE_PROPERTIES)+sizeof(KERNEL_LOGGER_NAME);
    ic.cswitchProperties = (EVENT_TRACE_PROPERTIES*)malloc(propertySize);
    memset(ic.cswitchProperties, 0, propertySize);
    // Fill the properties
    ic.cswitchProperties->EnableFlags      = EVENT_TRACE_FLAG_CSWITCH;   // That is what we want, context switches
    ic.cswitchProperties->LogFileMode      = EVENT_TRACE_REAL_TIME_MODE; // No file, retrieval through a callback
    ic.cswitchProperties->Wnode.BufferSize = propertySize;
    ic.cswitchProperties->Wnode.Flags      = WNODE_FLAG_TRACED_GUID;
    ic.cswitchProperties->Wnode.Guid       = SystemTraceControlGuid;
    ic.cswitchProperties->BufferSize       = 8; // In KB. Buffers are flushed when full, so small to be flushed often (else it is 1 Hz)
    ic.cswitchProperties->MinimumBuffers   = 1*PL_MAX_THREAD_QTY;
    ic.cswitchProperties->MaximumBuffers   = 4*PL_MAX_THREAD_QTY;
    ic.cswitchProperties->Wnode.ClientContext = 1;                 // 1 means rdtsc timer
    ic.cswitchProperties->LoggerNameOffset    = sizeof(EVENT_TRACE_PROPERTIES);
    memcpy(((char*)ic.cswitchProperties)+sizeof(EVENT_TRACE_PROPERTIES), KERNEL_LOGGER_NAME, sizeof(KERNEL_LOGGER_NAME));
    ic.cswitchPollEnabled = true; // Let's be optimistic

    // Stop the previous tracing (persistent across processes...) if not stopped properly by the last process using it.
    //  Of course, it modifies the property structure when the call really stops the previous tracing, hence the save/restore
    EVENT_TRACE_PROPERTIES save = *ic.cswitchProperties;
    ControlTrace(0, KERNEL_LOGGER_NAME, ic.cswitchProperties, EVENT_TRACE_CONTROL_STOP);
    *ic.cswitchProperties = save;

    // Start the tracing
    // Note: we will fail here if the executable is not run as administrator (not enough privileges)
    if(ic.cswitchPollEnabled && StartTrace(&ic.cswitchSessionHandle, KERNEL_LOGGER_NAME, ic.cswitchProperties)!=ERROR_SUCCESS) {
        ic.cswitchPollEnabled = false;
    }

    // Configure the logging. Indeed, tracing just activates the event collection so logging is required to retrieve the events...
    EVENT_TRACE_LOGFILE LogFile = {0};
#if defined(_UNICODE) || defined(UNICODE)
    LogFile.LoggerName          = (LPWSTR)KERNEL_LOGGER_NAME;
#else
    LogFile.LoggerName          = (LPSTR)KERNEL_LOGGER_NAME;
#endif
    LogFile.ProcessTraceMode    = PROCESS_TRACE_MODE_REAL_TIME | PROCESS_TRACE_MODE_EVENT_RECORD | PROCESS_TRACE_MODE_RAW_TIMESTAMP;
    LogFile.EventRecordCallback = plPriv::eventRecordCallback;

    if(ic.cswitchPollEnabled) {
        ic.cswitchConsumerHandle = OpenTrace(&LogFile);
        if(ic.cswitchConsumerHandle==(TRACEHANDLE)INVALID_HANDLE_VALUE) {
            CloseTrace(ic.cswitchSessionHandle);
            ic.cswitchPollEnabled = false;
        }
    }

    if(ic.cswitchPollEnabled) { // Successful activation
        ic.cswitchTraceLoggerThread = new std::thread(plPriv::collectCtxSwitch);
    } else free(ic.cswitchProperties); // Fail to activate the context switch logging
#endif // if defined(_WIN32) && PL_NOEVENT==0 && PL_IMPL_CONTEXT_SWITCH==1


    // Build the data exchange header
    int appNameLength   = (int)strlen(appName)+1;
    int buildNameLength = (buildName && buildName[0]!=0)? (int)strlen(buildName)+1 : 0;
#ifndef PL_PRIV_IMPL_LANGUAGE
#define PL_PRIV_IMPL_LANGUAGE "C++"
#endif
    int langNameLength  = strlen(PL_PRIV_IMPL_LANGUAGE)? (int)strlen(PL_PRIV_IMPL_LANGUAGE)+1 : 0;
    int tlvTotalSize    = 6/*Protocol TLV*/ + 20/*Time unit and origin*/ + 4+appNameLength/* Application name TLV */;
    if(buildNameLength) tlvTotalSize += 4+buildNameLength; // Optional build name TLV
    if(langNameLength)  tlvTotalSize += 4+langNameLength;  // Optional language name TLV
#if PL_EXTERNAL_STRINGS==1
    tlvTotalSize += 4;
#endif
#if PL_SHORT_STRING_HASH==1
    tlvTotalSize += 4;
#endif
#if PL_NOCONTROL==1
    tlvTotalSize += 4;
#endif
#if PL_SHORT_DATE==1
    tlvTotalSize += 12;
#endif
#if PL_COMPACT_MODEL==1
    tlvTotalSize += 4;
#endif
#if PL_IMPL_AUTO_INSTRUMENT==1
    ic.hasAutoInstrument = true;
#endif
    if(ic.hasAutoInstrument) {
        tlvTotalSize += 4;
    }
    if(ic.cswitchPollEnabled) {
        tlvTotalSize += 4;
    }
#if PL_HASH_SALT!=0
    tlvTotalSize += 8;
#endif
    int      headerSize = 16+tlvTotalSize;
    uint8_t* header = (uint8_t*)alloca(headerSize*sizeof(uint8_t));
    memset(header, 0, headerSize*sizeof(uint8_t)); // For the padding
    // Write the magic string to discriminate the connection type
    for(int i=0; i<8; ++i) header[i] = ("PL-MAGIC"[i]);
    // Write the endianess detection (provision, as little endian is supposed at the moment)
    *(uint32_t*)&header[8]  = 0x12345678;
    // Write the size of TLV block
    header[12] = (tlvTotalSize>>24)&0xFF;
    header[13] = (tlvTotalSize>>16)&0xFF;
    header[14] = (tlvTotalSize>> 8)&0xFF;
    header[15] = (tlvTotalSize    )&0xFF;
    // Write TLVs in big endian, T=2 bytes L=2 bytes
    int offset = 16;
    // TLV Protocol
    header[offset+0] = PL_TLV_PROTOCOL>>8; header[offset+1] = PL_TLV_PROTOCOL&0xFF;
    header[offset+2] = 0; header[offset+3] = 2; // 2 bytes payload
    header[offset+4] = PALANTEER_CLIENT_PROTOCOL_VERSION>>8; header[offset+5] = PALANTEER_CLIENT_PROTOCOL_VERSION&0xFF;
    offset += 6;
    // TLV Clock info
    header[offset+0] = PL_TLV_CLOCK_INFO>>8; header[offset+1] = PL_TLV_CLOCK_INFO&0xFF;
    header[offset+2] = 0; header[offset+3] = 16; // 16 bytes payload
    uint64_t clockOriginTick = PL_GET_CLOCK_TICK_FUNC();
    header[offset+ 4] = (clockOriginTick>>56)&0xFF; header[offset+ 5] = (clockOriginTick>>48)&0xFF;
    header[offset+ 6] = (clockOriginTick>>40)&0xFF; header[offset+ 7] = (clockOriginTick>>32)&0xFF;
    header[offset+ 8] = (clockOriginTick>>24)&0xFF; header[offset+ 9] = (clockOriginTick>>16)&0xFF;
    header[offset+10] = (clockOriginTick>> 8)&0xFF; header[offset+11] = (clockOriginTick    )&0xFF;
    char* tmp1 = (char*)&ic.tickToNs; uint64_t tmp = *(uint64_t*)tmp1;      // Avoids warning about strict aliasing
    header[offset+12] = (tmp>>56)&0xFF; header[offset+13] = (tmp>>48)&0xFF; // Standard IEEE 754 format, big endian
    header[offset+14] = (tmp>>40)&0xFF; header[offset+15] = (tmp>>32)&0xFF;
    header[offset+16] = (tmp>>24)&0xFF; header[offset+17] = (tmp>>16)&0xFF;
    header[offset+18] = (tmp>> 8)&0xFF; header[offset+19] = (tmp    )&0xFF;
    offset += 20;
    // TLV App name
    header[offset+0] = PL_TLV_APP_NAME>>8; header[offset+1] = PL_TLV_APP_NAME&0xFF;
    header[offset+2] = (appNameLength>>8)&0xFF; header[offset+3] = appNameLength&0xFF;
    memcpy(&header[offset+4], appName, appNameLength);
    offset += 4+appNameLength;
    // TLV build name
    if(buildNameLength>0) {
        header[offset+0] = PL_TLV_HAS_BUILD_NAME>>8; header[offset+1] = PL_TLV_HAS_BUILD_NAME&0xFF;
        header[offset+2] = (buildNameLength>>8)&0xFF; header[offset+3] = buildNameLength&0xFF;
        memcpy(&header[offset+4], buildName, buildNameLength);
        offset += 4+buildNameLength;
    }
    // TLV language name
    if(langNameLength>0) {
        header[offset+0] = PL_TLV_HAS_LANG_NAME>>8; header[offset+1] = PL_TLV_HAS_LANG_NAME&0xFF;
        header[offset+2] = (langNameLength>>8)&0xFF; header[offset+3] = langNameLength&0xFF;
        memcpy(&header[offset+4], PL_PRIV_IMPL_LANGUAGE, langNameLength);
        offset += 4+langNameLength;
    }

#define ADD_TLV_FLAG(flag)                                              \
    header[offset+0] = (flag)>>8; header[offset+1] = (flag)&0xFF;       \
    header[offset+2] = 0; header[offset+3] = 0; /* 0 bytes payload */   \
    offset += 4
#if PL_EXTERNAL_STRINGS==1
    ADD_TLV_FLAG(PL_TLV_HAS_EXTERNAL_STRING);
#endif
#if PL_SHORT_STRING_HASH==1
    ADD_TLV_FLAG(PL_TLV_HAS_SHORT_STRING_HASH);
#endif
#if PL_NOCONTROL==1
    ADD_TLV_FLAG(PL_TLV_HAS_NO_CONTROL);
#endif
#if PL_COMPACT_MODEL==1
    ADD_TLV_FLAG(PL_TLV_HAS_COMPACT_MODEL);
#endif
    if(ic.hasAutoInstrument) {
        ADD_TLV_FLAG(PL_TLV_HAS_AUTO_INSTRUMENT);
    }
    if(ic.cswitchPollEnabled) {
        ADD_TLV_FLAG(PL_TLV_HAS_CSWITCH_INFO);
    }
    plPriv::globalCtx.originNs = PL_GET_SYSTEM_CLOCK_NS();  // The global system date
#if PL_SHORT_DATE==1
    header[offset+ 0] = PL_TLV_HAS_SHORT_DATE>>8; header[offset+1] = PL_TLV_HAS_SHORT_DATE&0xFF;
    header[offset+ 2] = 0; header[offset+3] = 8; // 8 bytes payload
    tmp = plPriv::globalCtx.originNs;
    header[offset+ 4] = (tmp>>56)&0xFF; header[offset+ 5] = (tmp>>48)&0xFF;
    header[offset+ 6] = (tmp>>40)&0xFF; header[offset+ 7] = (tmp>>32)&0xFF;
    header[offset+ 8] = (tmp>>24)&0xFF; header[offset+ 9] = (tmp>>16)&0xFF;
    header[offset+10] = (tmp>> 8)&0xFF; header[offset+11] = (tmp    )&0xFF;
    ic.lastWrapCheckDateTick = (uint32_t)clockOriginTick;  // The associated event timestamp in tick
    ic.lastDateWrapQty       = 0;
    offset += 12;
#endif
#if PL_HASH_SALT!=0
    header[offset+0] = PL_TLV_HAS_HASH_SALT>>8; header[offset+1] = PL_TLV_HAS_HASH_SALT&0xFF;
    header[offset+2] = 0; header[offset+3] = 4; // 4 bytes payload
    header[offset+4] = (PL_HASH_SALT>>24)&0xFF; header[offset+5] = (PL_HASH_SALT>>16)&0xFF;
    header[offset+5] = (PL_HASH_SALT>> 8)&0xFF; header[offset+7] = (PL_HASH_SALT    )&0xFF;
    offset += 8;
#endif
    plAssert(offset==headerSize);

    // Write/send the built header
    bool headerSendingStatus = plPriv::palComSend(&header[0], headerSize);
    plAssert(headerSendingStatus, "Unable to send the session header");

    {
        // Create the transmission thread and wait for its readiness
        plAssert(PL_IMPL_STRING_BUFFER_BYTE_QTY>=128, "A minimum buffer size is required", PL_IMPL_STRING_BUFFER_BYTE_QTY);
        ic.txIsStarted = false;
        ic.threadServerTx = new std::thread(plPriv::transmitToServer);
        std::unique_lock<std::mutex> lk(ic.threadInitMx);
        ic.threadInitCv.wait(lk, [&] { return ic.txIsStarted; });
    }

#endif // if PL_NOCONTROL==0 || PL_NOEVENT==0
}


void
plStopAndUninit(void)
{
    auto& ic = plPriv::implCtx; PL_UNUSED(ic);

    // Unregister signals
#if PL_IMPL_CATCH_SIGNALS==1
    if(ic.signalHandlersSaved) {
        std::signal(SIGABRT, ic.signalsOldHandlers[0]);
        std::signal(SIGFPE,  ic.signalsOldHandlers[1]);
        std::signal(SIGILL,  ic.signalsOldHandlers[2]);
        std::signal(SIGSEGV, ic.signalsOldHandlers[3]);
        std::signal(SIGINT,  ic.signalsOldHandlers[4]);
        std::signal(SIGTERM, ic.signalsOldHandlers[5]);
#if defined(__unix__)
        std::signal(SIGPIPE, ic.signalsOldHandlers[6]);
#endif
#if defined(_WIN32)
        RemoveVectoredExceptionHandler(ic.exceptionHandler);
#endif // if defined(_WIN32)
    }
#endif // if PL_IMPL_CATCH_SIGNALS==1

#if PL_NOCONTROL==0 || PL_NOEVENT==0
    // Stop the data collection thread
    plPriv::globalCtx.enabled = false;
    {
        // Notify end of collection thread and wake it up
        std::lock_guard<std::mutex> lk(ic.txThreadSyncMx);
        ic.threadServerFlagStop.store(1);
        ic.txThreadSyncCv.notify_one();
    }
    if(ic.doNotUninit) {
        // Wait for the TX thread to send the last data sending (unless it is the crashing thread)
        if(ic.threadServerTx && ic.threadServerTx->joinable() && (int)PL_GET_SYS_THREAD_ID()!=ic.txThreadId) ic.threadServerTx->join();
        // No cleaning, so stop here
        return;
    }
    if(ic.threadServerTx && ic.threadServerTx->joinable()) ic.threadServerTx->join();
    if(ic.threadServerRx && ic.threadServerRx->joinable()) ic.threadServerRx->join();
#if defined(_WIN32) && PL_NOEVENT==0 && PL_IMPL_CONTEXT_SWITCH==1
    if(ic.cswitchTraceLoggerThread && ic.cswitchTraceLoggerThread->joinable()) {
        ic.cswitchTraceLoggerThread->join();
    }
#endif
    plPriv::globalCtx.collectEnabled = false;
    delete ic.threadServerTx; ic.threadServerTx = 0;
    delete ic.threadServerRx; ic.threadServerRx = 0;
#if defined(_WIN32) && PL_NOEVENT==0 && PL_IMPL_CONTEXT_SWITCH==1
    delete ic.cswitchTraceLoggerThread; ic.cswitchTraceLoggerThread = 0;
#endif

    // Restore the initial global state
    ic.threadServerFlagStop.store(0);
    delete[] ic.allocCollectBuffer; ic.allocCollectBuffer = 0;
    delete[] ic.sendBuffer; ic.sendBuffer = 0;
    plPriv::globalCtx.collectBuffers[0] = 0;
    plPriv::globalCtx.collectBuffers[1] = 0;
    plPriv::globalCtx.bankAndIndex.store(0);
    ic.prevBankAndIndex = 1UL<<31;
    ic.lkupStringToIndex.clear();
    ic.strBuffer.clear();
    ic.stringUniqueId = 0;
    ic.rxIsStarted = false;
    ic.txIsStarted = false;
    ic.frozenLastThreadBitmap = 0;
    ic.frozenThreadBitmap.store(0);
    ic.frozenThreadBitmapChange.store(0);
    ic.frozenThreadEnabled.store(0);

#endif // if PL_NOCONTROL==0 || PL_NOEVENT==0
}


plStats
plGetStats(void) { return plPriv::implCtx.stats; }


void
plDeclareVirtualThread(uint32_t externalVirtualThreadId, const char* format, ...)
{
    PL_UNUSED(externalVirtualThreadId); PL_UNUSED(format);
#if PL_NOEVENT==0
#if PL_VIRTUAL_THREADS==1

    // Build the name
    char name[PL_DYN_STRING_MAX_SIZE];
    va_list args;
    va_start(args, format);
    (void)vsnprintf(name, sizeof(name), format, args);
    va_end(args);

    // Ensure that the OS thread owns an internal Id (for context switches at least)
    plPriv::ThreadContext_t* tCtx = &plPriv::threadCtx;
    if(tCtx->id==0xFFFFFFFF) plPriv::getThreadId();

    // Get the internal thread Id
    plPriv::hashStr_t hash = (PL_FNV_HASH_OFFSET_^externalVirtualThreadId)*PL_FNV_HASH_PRIME_;
    uint32_t   newThreadId = 0;
    if(plPriv::implCtx.vThreadLkupExtToCtx.find(hash, newThreadId)) return; // Nothing to do as the vThread is already declared and only the first call matters

    // Create this new internal thread and virtual thread context
    newThreadId = plPriv::globalCtx.nextThreadId.fetch_add(1);
    plPriv::implCtx.vThreadLkupExtToCtx.insert(hash, newThreadId);

    if(newThreadId<PL_MAX_THREAD_QTY) {
        // Store the thread infos
        plPriv::ThreadInfo_t& ti = plPriv::globalCtx.threadInfos[newThreadId];
        ti.pid = 0xFFFFFFFF; // Prevent OS thread matching for context switches
        ti.nameHash = plPriv::hashString(name);
        int nameLength = (int)strlen(name);
        if(nameLength) memcpy(ti.name, name, nameLength+1);
        ti.name[PL_DYN_STRING_MAX_SIZE-1] = 0;

        // Log the thread declaration event with this new thread ID (not the standard call, as we do not want to override the real OS thread data)
        if(PL_IS_ENABLED_()) {
            uint32_t prevInternalThreadId = tCtx->id;
            tCtx->id = newThreadId;
            plPriv::eventLogRawDynName(PL_STRINGHASH(PL_BASEFILENAME), PL_EXTERNAL_STRINGS?0:PL_BASEFILENAME, name, 0, 0, PL_FLAG_TYPE_THREADNAME, 0);
            tCtx->id = prevInternalThreadId; // Restore previous thread Id
        }
    }

#else
    plAssert(0, "plDeclareVirtualThread is a specific API of the 'virtual threads' feature which is not enabled. Add -DPL_VIRTUAL_THREADS=1 in your compilation options.");
#endif // if PL_VIRTUAL_THREADS==1
#endif // if PL_NOEVENT==0
}


void
plDetachVirtualThread(bool isSuspended)
{
    PL_UNUSED(isSuspended);
#if PL_NOEVENT==0
#if PL_VIRTUAL_THREADS==1
    plPriv::ThreadContext_t* tCtx = &plPriv::threadCtx;

    // Ensure that the OS thread Id owns an internal Id (for context switches at least)
    if(tCtx->id==0xFFFFFFFF) plPriv::getThreadId();

    if(tCtx->id==tCtx->realId) return; // Nothing to do
    int vThreadId = tCtx->id;

    if(vThreadId<PL_MAX_THREAD_QTY) {
        // Save the suspend state
        plPriv::globalCtx.threadInfos[vThreadId].isSuspended = isSuspended;

        // Suspension is represented as an interrupting IRQ
        if(isSuspended && PL_IS_ENABLED_()) {
            plPriv::eventLogRaw(PL_STRINGHASH(""), PL_STRINGHASH("Suspended"), PL_EXTERNAL_STRINGS?0:"", "Suspended", 0, 0,
                                PL_FLAG_TYPE_SOFTIRQ | PL_FLAG_SCOPE_BEGIN, PL_GET_CLOCK_TICK_FUNC());
        }

        // The worker thread is no more used (logged on the virtual threadId)
        if(tCtx->realRscNameHash!=0 && PL_IS_ENABLED_()) {
            plPriv::eventLogRaw(PL_STRINGHASH(PL_BASEFILENAME), tCtx->realRscNameHash, PL_EXTERNAL_STRINGS?0:PL_BASEFILENAME, 0, 0, 0,
                                PL_FLAG_TYPE_LOCK_RELEASED, PL_GET_CLOCK_TICK_FUNC());
        }
    }

    // Switch to OS thread Id
    tCtx->id = tCtx->realId;

    if(vThreadId<PL_MAX_THREAD_QTY && PL_IS_ENABLED_()) {
        plPriv::ThreadInfo_t& ti = plPriv::globalCtx.threadInfos[vThreadId];
        if(ti.nameHash!=0 && ti.isBeginSent) {
            plPriv::eventLogRaw(PL_STRINGHASH(PL_BASEFILENAME), ti.nameHash, PL_EXTERNAL_STRINGS?0:PL_BASEFILENAME, 0, 0, 0,
                                PL_FLAG_SCOPE_END | PL_FLAG_TYPE_DATA_TIMESTAMP, PL_GET_CLOCK_TICK_FUNC());
            ti.isBeginSent = false;
        }
    }

#else  // if PL_VIRTUAL_THREADS==1
    plAssert(0, "plDetachVirtualThread is a specific API of the 'virtual thread' feature which is not enabled. Add -DPL_VIRTUAL_THREADS=1 in your compilation options.");
#endif // if PL_VIRTUAL_THREADS==1
#endif // if PL_NOEVENT==0
}


bool
plAttachVirtualThread(uint32_t externalVirtualThreadId)
{
    PL_UNUSED(externalVirtualThreadId);
    bool isNewVirtualThread = false;
#if PL_NOEVENT==0
#if PL_VIRTUAL_THREADS==1
    plPriv::ThreadContext_t* tCtx = &plPriv::threadCtx;

    // Ensure that the OS thread Id owns an internal Id (for context switches at least)
    if(tCtx->id==0xFFFFFFFF) plPriv::getThreadId();

    // Get the new virtual thread Id
    plPriv::hashStr_t hash = (PL_FNV_HASH_OFFSET_^externalVirtualThreadId)*PL_FNV_HASH_PRIME_;
    uint32_t vThreadId   = 0;
    if(!plPriv::implCtx.vThreadLkupExtToCtx.find(hash, vThreadId)) {
        // Create this new internal thread and virtual thread context
        vThreadId = plPriv::globalCtx.nextThreadId.fetch_add(1);
        plPriv::implCtx.vThreadLkupExtToCtx.insert(hash, vThreadId);
        isNewVirtualThread = true;
        if(vThreadId<PL_MAX_THREAD_QTY) {
            plPriv::globalCtx.threadInfos[vThreadId].pid = 0xFFFFFFFF; // Prevent OS thread matching for context switches
        }
    }
    if(tCtx->id==vThreadId) return isNewVirtualThread; // Nothing to do

    if(PL_IS_ENABLED_() && tCtx->realRscNameHash!=0 && tCtx->id!=tCtx->realId) {
        // The worker thread is switched to idle
        plPriv::eventLogRaw(PL_STRINGHASH(PL_BASEFILENAME), tCtx->realRscNameHash, PL_EXTERNAL_STRINGS?0:PL_BASEFILENAME, 0, 0, 0,
                            PL_FLAG_TYPE_LOCK_RELEASED, PL_GET_CLOCK_TICK_FUNC());
    }

    if(vThreadId<PL_MAX_THREAD_QTY && PL_IS_ENABLED_()) {
        plPriv::ThreadInfo_t& ti = plPriv::globalCtx.threadInfos[vThreadId];
        if(ti.nameHash!=0 && !ti.isBeginSent) {
            tCtx->id = tCtx->realId;
            plPriv::eventLogRaw(PL_STRINGHASH(PL_BASEFILENAME), ti.nameHash, PL_EXTERNAL_STRINGS?0:PL_BASEFILENAME, 0, 0, 0,
                                PL_FLAG_SCOPE_BEGIN | PL_FLAG_TYPE_DATA_TIMESTAMP, PL_GET_CLOCK_TICK_FUNC());
            ti.isBeginSent = true;
        }
    }

    // Overwrite the OS thread Id
    tCtx->id = vThreadId;

    // The worker thread is now dedicated to this virtual thread (logged on the worker threadId)
    if(vThreadId<PL_MAX_THREAD_QTY && plPriv::globalCtx.threadInfos[vThreadId].isSuspended && PL_IS_ENABLED_()) {
        plPriv::eventLogRaw(PL_STRINGHASH(""), PL_STRINGHASH("Suspended"), PL_EXTERNAL_STRINGS?0:"", "Suspended", 0, 0,
                            PL_FLAG_TYPE_SOFTIRQ | PL_FLAG_SCOPE_END, PL_GET_CLOCK_TICK_FUNC());
    }

    if(tCtx->realRscNameHash!=0 && PL_IS_ENABLED_()) {
        plPriv::eventLogRaw(PL_STRINGHASH(PL_BASEFILENAME), tCtx->realRscNameHash, PL_EXTERNAL_STRINGS?0:PL_BASEFILENAME, 0, 0, 0,
                            PL_FLAG_TYPE_LOCK_ACQUIRED, PL_GET_CLOCK_TICK_FUNC());
    }

#else  // if PL_VIRTUAL_THREADS==1
    plAssert(0, "plAttachVirtualThread is a specific API of the 'virtual thread' feature which is not enabled. Add -DPL_VIRTUAL_THREADS=1 in your compilation options.");
#endif // if PL_VIRTUAL_THREADS==1
#endif // if PL_NOEVENT==0
    return isNewVirtualThread;
}


void
plCrash(const char* message)
{
#if PL_NOCONTROL==0 || PL_NOEVENT==0
    // Do not log if the crash is located inside the Palanteer transmisison thread
    if(plPriv::implCtx.threadServerTx && (int)PL_GET_SYS_THREAD_ID()==plPriv::implCtx.txThreadId) {
        plPriv::globalCtx.enabled = false;
        plPriv::globalCtx.collectEnabled = false;
    }
#endif

    // Log and display the crash message
    plLogError("CRASH", "%s", message);
#if PL_IMPL_CONSOLE_COLOR==1
    PL_IMPL_PRINT_STDERR("\033[91m", true, false);  // Red
#endif
    PL_IMPL_PRINT_STDERR(message, true, false);
#if PL_IMPL_CONSOLE_COLOR==1
    PL_IMPL_PRINT_STDERR("\033[0m", true, false); // Standard
#endif
    PL_IMPL_PRINT_STDERR("\n", true, false);

    // Log and display the call stack
    PL_IMPL_STACKTRACE_FUNC();
    PL_IMPL_PRINT_STDERR("\n", true, true); // End of full crash display

    // Properly stop any recording, but without cleaning
    plPriv::implCtx.doNotUninit = true;
    plStopAndUninit();

    // Stop the process
    PL_IMPL_CRASH_EXIT_FUNC();
}



#endif //USE_PL==1

#undef PL_IMPLEMENTATION // 1

