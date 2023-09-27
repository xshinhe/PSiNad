#include "nadtcfer.h"

#include "../utils/definitions.h"

NAD_TCFer::NAD_TCFer(const int &itype,            // tcf_type
                     const int &nspec,            // get specification
                     const int &nsamp,            // get nsamp
                     const int &N, const int &F,  // problem size
                     const int &occ0              // help to get lentcf
) {
    tcf_type = itype;

    try {  // allocate bool matrix
        ALLOCATE_PTR_TO_VECTOR(tcf_0_bool, F * F);
        ALLOCATE_PTR_TO_VECTOR(tcf_t_bool, F * F);
        ALLOCATE_PTR_TO_VECTOR(tcf_0t_bool, F * F * F * F);
        ALLOCATE_PTR_TO_VECTOR(tcf_0_val, F * F);
        ALLOCATE_PTR_TO_VECTOR(tcf_t_val, F * F);
        ALLOCATE_PTR_TO_VECTOR(tcf_0t_val, F * F * F * F);
    } catch (std::runtime_error &e) { LOG(FATAL) << e.what(); }

    switch (tcf_type) {
        case nad_tcf::pop: {
            for (int i = 0, idx = 0; i < F; ++i)
                for (int j = 0; j < F; ++j, ++idx)
                    tcf_0_bool[idx] = (i == occ0 && j == occ0), tcf_t_bool[idx] = (i == j);
            lentcf = F;
            break;
        }
        case nad_tcf::rho: {
            for (int i = 0, idx = 0; i < F; ++i)
                for (int j = 0; j < F; ++j, ++idx) tcf_0_bool[idx] = (i == occ0 && j == occ0), tcf_t_bool[idx] = true;
            lentcf = F * F;
            break;
        }
        case nad_tcf::all: {
            for (int i = 0, idx = 0; i < F; ++i)
                for (int j = 0; j < F; ++j, ++idx) tcf_0_bool[idx] = true, tcf_t_bool[idx] = true;
            lentcf = F * F * F * F;
            break;
        }
        case nad_tcf::tcf:
        case nad_tcf::redtcf: {
            LOG(WARNING);
            try {  // read bool matrix from file
                int num1 = 0, num2 = 0;
                double tmp;
                std::ifstream ifs("tcf0");
                if (!ifs) LOG(FATAL) << "tcf file (0) cannot open!";

                for (int i = 0; i < F * F; ++i) {
                    if (ifs >> tmp) {
                        tcf_0_bool[i] = !(tmp == 0.0);
                        if (tcf_0_bool[i]) {
                            tcf_0_val[i] = tmp;
                            num1++;
                        }
                    }
                }
                ifs.close();

                ifs.open("tcft");
                if (!ifs) LOG(FATAL) << "tcf file (t) cannot open!";
                for (int i = 0; i < F * F; ++i) {
                    if (ifs >> tmp) {
                        tcf_t_bool[i] = !(tmp == 0.0);
                        if (tcf_t_bool[i]) {
                            tcf_t_val[i] = tmp;
                            num2++;
                        }
                    }
                }
                ifs.close();
                lentcf = num1 * num2;
            } catch (std::runtime_error &e) { LOG(FATAL) << e.what(); }
            break;
        }
        default:
            LOG(FATAL);
    }
    for (int i0 = 0, idx = 0; i0 < F * F; ++i0) {
        for (int it = 0; it < F * F; ++it, ++idx) {
            tcf_0t_bool[idx] = (tcf_0_bool[i0] && tcf_t_bool[it]);
            tcf_0t_val[idx]  = tcf_0_val[i0] * tcf_t_val[it];
        }
    }

    tcf_reduced = (tcf_type == nad_tcf::redtcf);
    if (tcf_reduced) lentcf = 1;

    CHECK_GT(lentcf, 0);

    this->nspec = nspec;
    this->lsamp = nspec * lentcf;
    this->nsamp = nsamp;
    this->ncoll = nsamp * lsamp;

    std::stringstream sss1, sss2;
    sss1 << FMT(0) << "stat" << FMT(8) << "time";
    if (tcf_reduced) {
        sss2 << FMT(8) << "Re(tcf)" << FMT(8) << "Im(tcf)";
    } else {
        for (int i0 = 0, i01 = 0; i01 < F; ++i01)
            for (int i02 = 0; i02 < F; ++i02, ++i0)
                for (int it = 0, it1 = 0; it1 < F; ++it1)
                    for (int it2 = 0; it2 < F; ++it2, ++it)
                        if (tcf_0_bool[i0] && tcf_t_bool[it]) {
                            sss2 << FMT(8) << utils::concat("Re(", i01, ",", i02, "|", it1, ",", it2, ")");
                            sss2 << FMT(8) << utils::concat("Im(", i01, ",", i02, "|", it1, ",", it2, ")");
                        }
    }
    header = sss1.str();
    for (int i = 0; i < nspec; ++i) header += sss2.str();

    try {
        ALLOCATE_PTR_TO_VECTOR(val, lentcf);
        ALLOCATE_PTR_TO_VECTOR(coll, ncoll);
        ALLOCATE_PTR_TO_VECTOR(stat, nsamp);
    } catch (std::bad_alloc &e) { LOG(FATAL) << e.what(); }
    Clear();
}

int NAD_TCFer::Clear() {
    for (int i = 0; i < nsamp; ++i) stat[i] = 0;
    for (int i = 0; i < ncoll; ++i) coll[i] = phys::math::iz;
    return 0;
}

int NAD_TCFer::Amount(NAD_TCFer &iT) {
    for (int i = 0; i < nsamp; ++i) stat[i] += iT.stat[i];
    for (int i = 0; i < ncoll; ++i) coll[i] += iT.coll[i];
    return 0;
}

int NAD_TCFer::MPIAmount(NAD_TCFer &iT) {
    Clear();
    MPI_Reduce(iT.coll, coll, ncoll, MPI::DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(iT.stat, stat, nsamp, MPI::INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return 0;
}

int NAD_TCFer::Count(const int &isamp, const int &now_spec) {
    stat[isamp] = 1;
    for (int ispec = 0, icount = lsamp * isamp, idx = 0; ispec < nspec; ++ispec) {
        if (ispec == now_spec)
            for (int i = 0; i < lentcf; ++i) coll[icount++] = val[i];
        else
            for (int i = 0; i < lentcf; ++i) coll[icount++] = phys::math::iz;
    }
    return 0;
}

bool NAD_TCFer::ifrecord(const int &i0, const int &it) { return (tcf_0_bool[i0] && tcf_t_bool[it]); }

int NAD_TCFer::report(const std::string &name, const double &unit) {
    std::ofstream ofs(name);
    ofs << std::setiosflags(std::ios::right) << header << std::endl;
    for (int i = 0, idx = 0; i < nsamp; ++i) {
        ofs << FMT(0) << stat[i] << FMT(8) << i * unit << " ";
        if (stat[i] > 0)
            for (int j = 0; j < lsamp; ++j, ++idx)
                ofs << FMT(8) << REAL_OF(coll[idx]) / stat[i]   // real part
                    << FMT(8) << IMAG_OF(coll[idx]) / stat[i];  // imag part
        else
            for (int j = 0; j < lsamp; ++j, ++idx) ofs << FMT(8) << REAL_OF(coll[idx]) << FMT(8) << IMAG_OF(coll[idx]);
        ofs << std::endl;
    }
    ofs.close();
    return 0;
}
