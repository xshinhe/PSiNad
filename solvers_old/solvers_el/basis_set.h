#ifndef BASIS_SET_H
#define BASIS_SET_H

#include <map>
#include <string>
#include <vector>

static const std::map<char, int> angular_map = {
    {'S', 0},
    {'P', 1},
    {'D', 2},
    {'F', 3},
};

enum class cart_btype { _S = 1, _P = 3, _SP = 4, _D = 6, _F = 10 };  ///< cartesian basis type
static const std::map<std::string, cart_btype> cart_btype_map = {
    {"S", cart_btype::_S}, {"P", cart_btype::_P}, {"SP", cart_btype::_SP}, {"D", cart_btype::_D}, {"F", cart_btype::_F},
};

enum class sphere_btype { _S = 1, _P = 3, _SP = 4, _D = 5, _F = 7 };  ///< spheric basis type
static const std::map<std::string, sphere_btype> sphere_btype_map = {
    {"S", sphere_btype::_S}, {"P", sphere_btype::_P}, {"SP", sphere_btype::_SP},
    {"D", sphere_btype::_D}, {"F", sphere_btype::_F},
};

struct Qnum {
    int L, M, N;
};

static const std::map<cart_btype, std::vector<Qnum>> qnums_map = {
    {cart_btype::_S, {{0, 0, 0}}},
    {cart_btype::_P, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}},
    {cart_btype::_SP, {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}},
    {cart_btype::_D, {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}}},
    {cart_btype::_F,
     {{3, 0, 0}, {0, 3, 0}, {0, 0, 3}, {2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 0, 2}, {0, 2, 1}, {0, 1, 2}, {1, 1, 1}}},
};

class Atomic_BasisSet {  // handle different basis
   public:
    using Btype = cart_btype;

    int znum;
    std::string bname;
    std::vector<Btype> rbtype;
    std::vector<double> rdata;
    std::vector<std::string> rdata_angl;
    std::vector<int> rdata_posi;
    std::vector<int> rdata_posf;
    std::vector<int> rdata_nrow;
    std::vector<int> rdata_ncol;
    int nb;  // ncgto
    int ngto;

    double* coeffs;
    double* alphas;
    Qnum* qnums;
    int* nconts;

    Atomic_BasisSet(const int& znum, const std::string& bname) : znum{znum}, bname{bname} {
        std::string basis_dir = ".";
        std::ifstream fin(bname);
        std::string eachline;
        std::string tmpbtype, tmp;
        double tmpd;

        if (!fin) LOG(FATAL) << "No basis file here {" << bname << "}";

        int cnt = 1;
        while (getline(fin, eachline)) {
            if (eachline[0] == '!') continue;
            if (eachline[0] == '*') cnt++;
            if (cnt < znum) continue;
            fin >> tmp >> tmp;  // name & charge
            while (fin >> tmpbtype) {
                if (tmpbtype[0] == '*') break;

                // else, tmpbtype is in ["S", "SP", "P", "D", "F", ...]
                // LOG(WARNING) << tmpbtype;

                rdata_angl.push_back(tmpbtype);

                rbtype.push_back(cart_btype_map.at(tmpbtype));
                rdata_posi.push_back(rdata.size());

                int nrow;
                fin >> nrow >> tmp;
                getline(fin, eachline);  // remove last "\n"

                rdata_nrow.push_back(nrow);

                for (int i = 0; i < nrow; ++i) {
                    getline(fin, eachline);
                    std::stringstream ss(eachline);
                    // LOG(WARNING) << eachline;

                    int ncol = 0;
                    while (ss >> tmpd) {  // each columns
                        LOG(WARNING) << tmpd;
                        rdata.push_back(tmpd);
                        if (i == 0) ncol++;
                    }
                    if (i == 0) rdata_ncol.push_back(ncol);
                }
                rdata_posf.push_back(rdata.size());
            }
            break;
        }
        fin.close();

        nb = 0, ngto = 0;
        for (int i = 0; i < rbtype.size(); ++i) {
            nb += (int) rbtype[i];
            ngto += ((int) rbtype[i]) * rdata_nrow[i];
        }

        nconts = new int[nb];
        qnums  = new Qnum[nb];
        alphas = new double[ngto];
        coeffs = new double[ngto];

        for (int i = 0, idx = 0, jdx = 0; i < rbtype.size(); ++i) {
            Btype& irbtype = rbtype[i];
            switch (irbtype) {
                case Btype::_SP: {
                    const std::vector<Qnum>& qnum_type_ref = qnums_map.at(irbtype);
                    for (int iq = 0; iq < (int) Btype::_SP; ++iq) {
                        nconts[jdx]  = rdata_nrow[i];
                        qnums[jdx++] = qnum_type_ref[iq];
                        for (int j = rdata_posi[i]; j < rdata_posf[i];) {
                            alphas[idx] = rdata[j++];
                            if (iq != 0) j++;  // skip S
                            coeffs[idx++] = rdata[j++];
                            if (iq == 0) j++;  // skip P
                        }
                    }
                    break;
                }
                case Btype::_S:
                case Btype::_P:
                case Btype::_D:
                case Btype::_F: {
                    const std::vector<Qnum>& qnum_type_ref = qnums_map.at(irbtype);
                    for (int iq = 0; iq < (int) irbtype; ++iq) {
                        nconts[jdx]  = rdata_nrow[i];
                        qnums[jdx++] = qnum_type_ref[iq];
                        for (int j = rdata_posi[i]; j < rdata_posf[i];) {
                            alphas[idx]   = rdata[j++];
                            coeffs[idx++] = rdata[j++];
                        }
                    }
                    break;
                }
                default: {
                    throw std::runtime_error("unsupported type now");
                    exit(-1);
                }
            }
        }
    }

    /* used by vector.push_back() (copy constructor without noexcept) */
    Atomic_BasisSet(const Atomic_BasisSet& ABS)
        : znum{ABS.znum},
          bname{ABS.bname},
          rbtype{ABS.rbtype},
          rdata{ABS.rdata},
          rdata_angl{ABS.rdata_angl},
          rdata_posi{ABS.rdata_posi},
          rdata_posf{ABS.rdata_posf},
          rdata_nrow{ABS.rdata_nrow},
          rdata_ncol{ABS.rdata_ncol},
          nb{ABS.nb},
          ngto{ABS.ngto} {
        nconts = new int[nb];
        qnums  = new Qnum[nb];
        alphas = new double[ngto];
        coeffs = new double[ngto];

        for (int i = 0; i < nb; ++i) { nconts[i] = ABS.nconts[i], qnums[i] = ABS.qnums[i]; }
        for (int i = 0; i < ngto; ++i) { alphas[i] = ABS.alphas[i], coeffs[i] = ABS.coeffs[i]; }
    }

    virtual ~Atomic_BasisSet() {
        delete[] nconts;
        delete[] coeffs;
        delete[] alphas;
        delete[] qnums;
    }
};



#endif  // BASIS_SET_H
