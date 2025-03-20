#include "kids/VariableDescriptor.h"

#include <iostream>
#include <regex>

#include "kids/debug_utils.h"
#include "kids/Types.h"

// if you want to debug the meta-information for the VariableDescriptor,
// uncomment the following line
// #define LOCAL_DEBUG
#ifdef LOCAL_DEBUG
#include "kids/fmt.h"
#endif // LOCAL_DEBUG

namespace PROJECT_NS {

/**
 * @class VariableDescriptor
 * @brief Represents a variable descriptor with metadata for parsing and
 * processing.
 *
 * This class parses a token string to extract variable information such as
 * name, field, time, index, and type. It also sets default values and derives
 * additional properties based on the input.
 */
VariableDescriptor::VariableDescriptor(const std::string& token_string,
                                       const std::string& save,
                                       VariableDescriptorPolicy::_type vtype_in)
    : tokenString{token_string}, // The original token string
      save{save},                // The save string (purpose to be defined)
      vtype{vtype_in}
{ // Flag indicating if the variable is a tabular output

    /**
     * Regular expression pattern to parse the token string for a variable
     * An example is:
     *
     *      `name{field@time}<index>:type`
     *
     * Here name is required and others are optional
     * Pattern breakdown: [note '(...)' is for a capturing group, '(?:...)' is
     * for a non-capturing group]
     * ([^{<:]*)        Capturing group 1 -> march name with any characters
     * except '{', '<', ':' from beginning
     * (?:              begin of a non-capturing group A
     * \\{([^@<]*)      capturing group 2 -> match field with '{' followed by
     * any characters except '@' and '<'
     * (?:@(.*?))?      Optional non-capturing group B and capturing group 3 ->
     *                                       match time with '@' followed by any
     * characters (non-greedy)
     * \})?             '\}' should concide with group 2 otherwise group A is
     * optional
     * (?:<([^>]*)>)?   Optional non-capturing group C and capturing group 4 ->
     *                                       match index with any characters
     * between '<' and '>' except '>'
     * (?::([^:]*))?    Optional non-capturing group D and capturing group 5 ->
     *                                       match type with ':' followed by any
     * characters except ':'
     */
    const std::regex pattern(
        "([^{<:]*)(?:\\{([^@<]*)(?:@(.*?))?\\})?(?:<([^>]*)>)?(?::([^:]*))?");
    std::smatch match; // Object to store the results of the regex match

    // Attempt to match the token string against the regex pattern
    if (std::regex_match(tokenString, match, pattern))
    {
        // match[1] to match[5] correspond to the capture groups in the regex
        std::tie(name, field, time, index, type) =
            std::make_tuple(match[1].str(), // Name of the variable
                            match[2].str(), // Field (optional)
                            match[3].str(), // Time (optional)
                            match[4].str(), // Index (optional)
                            match[5].str()  // Type (optional)
            );
    }
    else
    {
        // If the token string does not match the pattern, throw an error
        throw kids_error(
            utils::concat("Cannot match variable pattern: ", token_string));
    }

    // Process the 'field' to determine the namespace
    if (field.empty() || field == "I")
    {
        field = "integrator"; // Default namespace if field is empty or "I"
                              // (Shortcut for "integrator" namespace)
    }
    if (field == "M")
    {
        field = "model"; // Shortcut for "model" namespace
    }
    if (field == "P")
    {
        field = "parameter"; // Shortcut for "parameter" namespace
    }
    if (field == "R")
    {
        field = "record"; // Shortcut for "record" namespace
    }
    // if TabularOutput & MetaOutput
    if (vtype != VariableDescriptorPolicy::Input) field = "record";

    // Process the 'index' if present
    // Currently, if index is empty, do nothing
    // Further processing can be added here if needed
    if (index.empty())
    {
        // No action needed
    }

    // Process the 'type' to determine the data type
    dataType = kids_void_type; // Default data type is void
    if (type == "R")
    {
        dataType = kids_real_type; // Real number type
    }
    if (type == "C")
    {
        dataType = kids_complex_type; // Complex number type
    }
    if (type == "I")
    {
        dataType = kids_int_type; // Integer number type
        // untested
        throw kids_error("TODO");
    }

    // Process the 'name' and determine if the tabular variable shoube be
    // degrade to MetaOutput
    if (vtype == VariableDescriptorPolicy::TabularOutput && !name.empty() &&
        name[0] == '*')
    {
        vtype = VariableDescriptorPolicy::MetaOutput;
        name = name.substr(1); // Remove the '*' from the name
    }

    // Process the keys based on the variable properties
    keyRaw = utils::concat(
        field, ".", name); // Raw key is a concatenation of field and name
    keyRec = keyRaw;       // Record key defaults to raw key
    keyRes0 = keyRaw;      // Result keys default to raw key
    keyRes1 = keyRaw;
    keyRes2 = keyRaw;

    if (vtype == VariableDescriptorPolicy::TabularOutput)
    {
        // If the variable is an output, set the record and result keys
        // accordingly
        keyRec = utils::concat("record.", name); // Record key for output
        keyRes0 = utils::concat("_.0.", field, ".",
                                name); // Result key for first result
        keyRes1 = utils::concat("_.1.", field, ".",
                                name); // Result key for second result
        keyRes2 = utils::concat("_.2.", field, ".",
                                name); // Result key for third result
    }
    else // @TODO: so {@} cannot be used for TabularOutput
    {
        // If the variable is not an output, set the record key based on time if
        // present
        if (time.empty())
        {
            // if time namespace if not specified, keyRec is just keyRaw, i.e.,
            // only by reference
            keyRec =
                keyRaw; // No time specified, record key is the same as raw key
        }
        else
        {
            // otherwise, keyRec will copy keyRaw from a different storage
            // and the copy only occurs at time = sampIndex
            keyRec = utils::concat("@.", time, ".", field, ".",
                                   name); // Include time in record key
        }
    }
}

/**
 * @brief Defines the variable within a given DataSet with specified properties.
 *
 * This method registers the VariableDescriptor in the provided DataSet,
 * handling both input and output variables appropriately. It sets up the data
 * pointers and shapes based on the variable's type and whether it is tabular.
 *
 * @param DS The DataSet in which to define the variable.
 * @param data_type The data type of the variable (e.g., real, complex).
 * @param cxxshape The shape of the variable's data in C++ convention (row-major
 * order).
 * @param totalFrameNumber The total number of frames for tabular data.
 *
 * @throws kids_error If there is a conflict or loss of the variable in the
 * DataSet.
 */
void VariableDescriptor::defineIn(std::shared_ptr<DataSet> DS,
                                  kids_dtype data_type,
                                  const std::vector<std::size_t>& cxxshape,
                                  std::size_t totalFrameNumber)
{
    // For the rule expression `A(B,C,D)`
    switch (vtype)
    {
        case VariableDescriptorPolicy::TabularOutput:
        case VariableDescriptorPolicy::MetaOutput: {
            // Output variables must not be defined multiple times in the
            // DataSet
            // uncomment for check & debug
            // if (DS->haskey(keyRec)) {
            //     throw kids_error(utils::concat("Conflict of key: ", keyRec));
            // }

            // Initialize the raw data pointer to null, because it has no
            // reference storage
            dataPointerRaw = nullptr;

            // Prepare the shape for stacked data (e.g., tabular data)
            // if we have n frames for recording a tensor with the shape <abcd>,
            // we should prepare <nabcd> memory to store the values.
            std::vector<std::size_t> cxxstackedshape;
            cxxstackedshape.push_back(
                totalFrameNumber); // Prepend the total number of frames
            for (const auto& dim : cxxshape)
            {
                cxxstackedshape.push_back(
                    dim); // Append the original dimensions
            }

            // Define the variable in the DataSet based on its data type
            switch (data_type)
            {
                case kids_real_type: {
                    // Define a real-valued variable with the given key and
                    // shape
                    DS->def_real(keyRec, cxxshape,
                                 utils::concat(name, " traced in 1 frame"));
                    if (vtype == VariableDescriptorPolicy::TabularOutput)
                    { // save nframes of a tensor
                        // If the variable is tabular, define additional keys
                        // for results
                        DS->def_real(keyRes0, cxxstackedshape,
                                     utils::concat(
                                         name, " collected n frame of 1 traj"));
                        DS->def_real(
                            keyRes1, cxxstackedshape,
                            utils::concat(
                                name,
                                " collected n frame of m traj but 1 mpi"));
                        DS->def_real(
                            keyRes2, cxxstackedshape,
                            utils::concat(name, " collected n frame of k mpi"));
                    }
                    break;
                }

                case kids_complex_type: {
                    // Define a complex-valued variable with the given key and
                    // shape
                    DS->def_complex(keyRec, cxxshape,
                                    utils::concat(name, " traced in 1 frame"));
                    if (vtype == VariableDescriptorPolicy::TabularOutput)
                    { // save nframes of a tensor
                        // If the variable is tabular, define additional keys
                        // for results
                        DS->def_complex(
                            keyRes0, cxxstackedshape,
                            utils::concat(name,
                                          " collected n frame of 1 traj"));
                        DS->def_complex(
                            keyRes1, cxxstackedshape,
                            utils::concat(
                                name,
                                " collected n frame of m traj but 1 mpi"));
                        DS->def_complex(
                            keyRes2, cxxstackedshape,
                            utils::concat(name, " collected n frame of k mpi"));
                    }
                    break;
                }

                default: {
                    // Handle other data types if necessary
                    throw kids_error(utils::concat(
                        "Unsupported data type for output variable '", //
                        name, "': ", enum_t_as_str(dataType)));
                }
            }

            // Obtain the data type, pointer, and shape from the DataSet for the
            // main key
            std::tie(dataType, dataPointerTrace, shape) = DS->obtain(keyRec);
            if (vtype == VariableDescriptorPolicy::TabularOutput)
            {
                // If the variable is tabular, obtain additional pointers for
                // the result keys
                std::tie(dataType, dataPointerRes0, stackedshape) =
                    DS->obtain(keyRes0);
                std::tie(dataType, dataPointerRes1, stackedshape) =
                    DS->obtain(keyRes1);
                std::tie(dataType, dataPointerRes2, stackedshape) =
                    DS->obtain(keyRes2);
            }
            else
            {
                dataPointerRes0 = nullptr;
                dataPointerRes1 = nullptr;
                dataPointerRes2 = nullptr;
            }
            break;
        }

        case VariableDescriptorPolicy::Input: {
            // Input variables must already exist in the DataSet!
            if (!DS->haskey(keyRaw))
            {
                throw kids_error(utils::concat("Loss of key: ", keyRaw));
            }

            // for variables not of isOutput, keyRaw [->dataPointerRaw] &
            // keyRec[->dataPointerTrace] refer to the same storage

            // Obtain the data type, pointer, and shape from the DataSet for the
            // raw key
            std::tie(dataType, dataPointerRaw, shape) = DS->obtain(keyRaw);

            // Define the trace key based on the data type
            switch (dataType)
            {
                case kids_real_type: {
                    // Define a real-valued trace key with the same shape as the
                    // raw data
                    dataPointerTrace = (void*)DS->def_real(
                        keyRec, *shape, " traced in 1 frame");
                    break;
                }

                case kids_complex_type: {
                    // Define a complex-valued trace key with the same shape as
                    // the raw data
                    dataPointerTrace = (void*)DS->def_complex(
                        keyRec, *shape, " traced in 1 frame");
                    break;
                }

                default: {
                    // Handle other data types if necessary
                    throw kids_error(utils::concat(
                        "Unsupported data type for input variable '", //
                        name, "': ", enum_t_as_str(dataType)));
                }
            }

            // Initialize the stacked shape to null for inside variables
            dataPointerRes0 = nullptr;
            dataPointerRes1 = nullptr;
            dataPointerRes2 = nullptr;
            stackedshape = nullptr;

            break;
        }
    }
#ifdef LOCAL_DEBUG
    std::cout << LOC() << "VariableDescriptor with Info:\n"
              << ".tokenString = " << tokenString << "\n"          //
              << ".name = " << name << "\n"                        //
              << ".field = " << field << "\n"                      //
              << ".index = " << index << "\n"                      //
              << ".type = " << type << "\n"                        //
              << ".time = " << time << "\n"                        //
              << ".vtype = " << vtype << "\n"                      //
              << ".dataType = " << enum_t_as_str(dataType) << "\n" //
              << ".keyRaw = " << keyRaw << "\n"                    //
              << ".keyRec = " << keyRec << "\n"                    //
              << ".keyRes0 = " << keyRes0 << "\n"                  //
              << ".keyRes1 = " << keyRes1 << "\n"                  //
              << ".keyRes2 = " << keyRes2 << "\n";                 //
#endif // LOCAL_DEBUG
#ifdef LOCAL_DEBUG
    std::cout << LOC() << "VariableDescriptor with Data Storage:\n"
              << ".tokenString = " << tokenString << "\n"           //
              << ".dataType = " << enum_t_as_str(dataType) << "\n"  //
              << ".dataPointerRaw = " << dataPointerRaw << "\n"     //
              << ".dataPointerTrace = " << dataPointerTrace << "\n" //
              << ".dataPointerRes0 = " << dataPointerRes0 << "\n"   //
              << ".dataPointerRes1 = " << dataPointerRes1 << "\n"   //
              << ".dataPointerRes2 = " << dataPointerRes2 << "\n"   //
              << ".shape = " << ((shape) ? shape->to_string() : "") << "\n";
#endif // LOCAL_DEBUG
}

/**
 * @brief Checks and updates the trace data based on the current sample index if
 * time namespace exists.
 *
 * This method is responsible for copying data from the raw data pointer to the
 * trace data pointer if the current sample index is equal to the variable's
 * time index. It handles both real and complex data types.
 *
 * @param sampleIndex The current sample index to check against the variable's
 * time index.
 *
 * @throws kids_error If the data type is unsupported or if memory access is
 * invalid.
 */
void VariableDescriptor::checkTrace(int sampleIndex)
{
    // If the variable has no time index or is an output, there is no need to
    // check the trace because an inside variable without time namespace is
    // traced by reference (let keyRec = keyRaw) and outside output variable is
    // updated with only the rule.
    if (time.empty() || vtype != VariableDescriptorPolicy::Input) return;

    // Attempt to convert the time string to an integer
    int timeIndex;
    try
    {
        timeIndex = std::stoi(time);
    }
    catch (const std::invalid_argument& e)
    {
        throw kids_error(utils::concat("Invalid time index for variable '",
                                       name, "': ", time));
    }
    catch (const std::out_of_range& e)
    {
        throw kids_error(utils::concat("Time index out of range for variable '",
                                       name, "': ", time));
    }

    // If the current sample index matches the variable's time index, proceed to
    // copy data
    if (timeIndex == sampleIndex)
    {
        switch (dataType)
        {
            case kids_real_type: {
                // Cast the raw data pointer to a pointer of type kids_real
                kids_real* fromdata = static_cast<kids_real*>(dataPointerRaw);
                // Cast the trace data pointer to a pointer of type kids_real
                kids_real* todata = static_cast<kids_real*>(dataPointerTrace);

                // Check if shape is valid
                if (!shape)
                {
                    throw kids_error(utils::concat(
                        "Shape is null for variable '", name, "'"));
                }

                // Copy data from fromdata to todata
                for (std::size_t i = 0; i < shape->size(); ++i)
                {
                    todata[i] = fromdata[i];
                }
                break;
            }

            case kids_complex_type: {
                // Cast the raw data pointer to a pointer of type kids_complex
                kids_complex* fromdata =
                    static_cast<kids_complex*>(dataPointerRaw);
                // Cast the trace data pointer to a pointer of type kids_complex
                kids_complex* todata =
                    static_cast<kids_complex*>(dataPointerTrace);

                // Check if shape is valid
                if (!shape)
                {
                    throw kids_error(utils::concat(
                        "Shape is null for variable '", name, "'"));
                }

                // Copy data from fromdata to todata
                for (std::size_t i = 0; i < shape->size(); ++i)
                {
                    todata[i] = fromdata[i];
                }
                break;
            }

            default:
                // If the data type is unsupported, throw an error
                throw kids_error(utils::concat(
                    "Unsupported data type for variable '", name, "'"));
        }
    }
}

}; // namespace PROJECT_NS