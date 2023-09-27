#include "solver.h"

#include "../utils/definitions.h"
#include "../utils/hdf5_utils.h"

#define OPENDF_MPI_PARALLEL

Solver::Solver(const Param& iparm, Model* pM) : Model(iparm) {
    pModel = pM;
    if (pModel->iou.time != 1 || pModel->iou.leng != 1 || pModel->iou.time != 1 || pModel->iou.mass != 1 ||
        pModel->iou.ener != 1 || pModel->iou.temp != 1) {
        iou = pModel->iou;
    }
    para_flag = FLAGS_para_flag;
}

Solver::~Solver() { pModel = nullptr; }

int Solver::run() {
    para_type = ParaPolicy::_dict.at(para_flag);
    switch (para_type) {
        case ParaPolicy::no_para: {
            if (FLAGS_r == "") {
                pContext_in = new File("in.h5", File::ReadWrite | File::Create | File::Truncate);
                rand_rng = SeededEngine(nullptr);
            } else {
                if (mpi_rank == 0) utils::copyfile_from_to(FLAGS_r, "in.h5");
                try {
                    pContext_in = new File("in.h5", File::ReadWrite);
                } catch (HighFive::Exception& e) { LOG(FATAL); }
                if (FLAGS_read_seed) {
                    LOAD(pContext_in, rng_seeds, rng_t::state_size);
                    rand_rng = SeededEngine(rng_seeds);
                } else {
                    rand_rng = SeededEngine(nullptr);
                }
            }
            pContext_out = new File("out.h5", File::ReadWrite | File::Create | File::Truncate);
            run_impl();

            // save random seeds
            DUMP(pContext_out, rng_seeds, rng_t::state_size);
            
            pContext_out->flush();
            pContext_in->flush();
            delete pContext_out; delete pContext_in;
            break;
        }
        case ParaPolicy::calc_para: {
            mpi_utils_init();
            if (FLAGS_r == "") {
                pContext_in = new File("in.h5", File::ReadWrite | File::Create | File::Truncate,
                                    MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));

                rand_rng = SeededEngine(nullptr);
            } else {
                if (mpi_rank == 0) utils::copyfile_from_to(FLAGS_r, "rin.h5");
                try {
                    pContext_in = new File("in.h5", File::ReadWrite, MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
                } catch (HighFive::Exception& e) { LOG(FATAL); }

                if (FLAGS_read_seed) {
                    LOAD(pContext_in, rng_seeds, rng_t::state_size);
                    rand_rng = SeededEngine(rng_seeds);
                } else {
                    rand_rng = SeededEngine(nullptr);
                }
            }
            pContext_out = new File("out.h5", File::ReadWrite | File::Create | File::Truncate,
                                    MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
            run_impl();


            // save random seeds
            DUMP(pContext_out, rng_seeds, rng_t::state_size);
            
            pContext_out->flush();
            pContext_in->flush();
            delete pContext_out; delete pContext_in;

            mpi_utils_final();
            break;
        }
        case ParaPolicy::traj_para: {
            mpi_utils_init();
            if (FLAGS_r == "") {
                pContext_in = new File("in.h5", File::ReadWrite | File::Create | File::Truncate,
                                    MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
                rand_rng = SeededEngine(nullptr);
            } else {
                if (mpi_rank == 0) utils::copyfile_from_to(FLAGS_r, "in.h5");
                try {
                    pContext_in = new File("in.h5", File::ReadWrite, MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
                } catch (HighFive::Exception& e) { LOG(FATAL); }

                if (FLAGS_read_seed) {
                    LOAD(pContext_in, rng_seeds, rng_t::state_size);
                    rand_rng = SeededEngine(rng_seeds);
                } else {
                    rand_rng = SeededEngine(nullptr);
                }
            }
            pContext_out = new File("out.h5", File::ReadWrite | File::Create | File::Truncate,
                                    MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));

            run_parallel();

            // save random seeds
            DUMP(pContext_out, rng_seeds, rng_t::state_size);
            
            pContext_out->flush();
            pContext_in->flush();
            delete pContext_out; delete pContext_in;

            mpi_utils_final();
            break;
        }
    }
    return 0;
}

int Solver::run_impl() {
    LOG(FATAL);
    return 0;
}

int Solver::run_parallel() {
    LOG(FATAL);
    return 0;
}
