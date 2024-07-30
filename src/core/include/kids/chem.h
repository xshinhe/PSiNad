/**@file        chem.h
 * @brief       provide databases in chemistry
 * @details     all utils are provided under the namespace chem.
 *
 * @author      Xin He
 * @date        2024-04
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of Xin's PhD research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @warning    Do not include this file to any header. You'd better include it only
 *  in source files!
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-03-29  <td> initial version.
 * <tr><td> 2024-06-20  <td> updated version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef CHEM_H
#define CHEM_H
#include <algorithm>
#include <array>
#include <cstring>
#include <string>

namespace chem {

struct ElemInfo {
    ElemInfo(const std::string& label, int Z, int A, double mass_in_amu)  //
        : label{label}, Z{Z}, A{A}, mass{mass_in_amu} {};
    std::string label;
    int         Z;  // proton number
    int         A;  // proton + neotron number; A should greater than Z, otherwise A=-1 means mixed value, A=-1 means
                    // max-abundace value
    double mass;
};

namespace elem {
// clang-format off
const ElemInfo NU("NU",  0, -1, 0.0          );
const ElemInfo H ("H" ,  1, -1, 1.007947     );
const ElemInfo HE("HE",  2, -1, 4.0026022    );
const ElemInfo LI("LI",  3, -1, 6.9412       );
const ElemInfo BE("BE",  4, -1, 9.0121823    );
const ElemInfo B ("B" ,  5, -1, 10.8117      );
const ElemInfo C ("C" ,  6, -1, 12.01078     );
const ElemInfo N ("N" ,  7, -1, 14.00672     );
const ElemInfo O ("O" ,  8, -1, 15.99943     );
const ElemInfo F ("F" ,  9, -1, 18.99840325  );
const ElemInfo NE("NE", 10, -1, 20.17976     );
const ElemInfo NA("NA", 11, -1, 22.989769282 );
const ElemInfo MG("MG", 12, -1, 24.30506     );
const ElemInfo AL("AL", 13, -1, 26.98153868  );
const ElemInfo SI("SI", 14, -1, 28.08553     );
const ElemInfo P ("P" , 15, -1, 30.9737622   );
const ElemInfo S ("S" , 16, -1, 32.0655      );
const ElemInfo CL("CL", 17, -1, 35.4532      );
const ElemInfo AR("AR", 18, -1, 39.9481      );
const ElemInfo K ("K" , 19, -1, 39.09831     );
const ElemInfo CA("CA", 20, -1, 40.0784      );
const ElemInfo SC("SC", 21, -1, 44.9559126   );
const ElemInfo TI("TI", 22, -1, 47.8671      );
const ElemInfo V ("V" , 23, -1, 50.94151     );
const ElemInfo CR("CR", 24, -1, 51.99616     );
const ElemInfo MN("MN", 25, -1, 54.9380455   );
const ElemInfo FE("FE", 26, -1, 55.8452      );
const ElemInfo CO("CO", 27, -1, 58.9331955   );
const ElemInfo NI("NI", 28, -1, 58.69342     );
const ElemInfo CU("CU", 29, -1, 63.5463      );
const ElemInfo ZN("ZN", 30, -1, 65.4094      );
const ElemInfo GA("GA", 31, -1, 69.7231      );
const ElemInfo GE("GE", 32, -1, 72.641       );
const ElemInfo AS("AS", 33, -1, 74.921602    );
const ElemInfo SE("SE", 34, -1, 78.963       );
const ElemInfo BR("BR", 35, -1, 79.9041      );
const ElemInfo KR("KR", 36, -1, 83.7982      );
const ElemInfo RB("RB", 37, -1, 85.46783     );
const ElemInfo SR("SR", 38, -1, 87.621       );
const ElemInfo Y ("Y" , 39, -1, 88.905852    );
const ElemInfo ZR("ZR", 40, -1, 91.2242      );
const ElemInfo NB("NB", 41, -1, 92.906382    );
const ElemInfo MO("MO", 42, -1, 95.942       );
const ElemInfo TC("TC", 43, -1, 98.0         );
const ElemInfo RU("RU", 44, -1, 101.072      );
const ElemInfo RH("RH", 45, -1, 102.905502   );
const ElemInfo PD("PD", 46, -1, 106.421      );
const ElemInfo AG("AG", 47, -1, 107.86822    );
const ElemInfo CD("CD", 48, -1, 112.4118     );
const ElemInfo IN("IN", 49, -1, 114.8183     );
const ElemInfo SN("SN", 50, -1, 118.7107     );
const ElemInfo SB("SB", 51, -1, 121.7601     );
const ElemInfo TE("TE", 52, -1, 127.603      );
const ElemInfo I ("I" , 53, -1, 126.904473   );
const ElemInfo XE("XE", 54, -1, 131.2936     );
const ElemInfo CS("CS", 55, -1, 132.90545192 );
const ElemInfo BA("BA", 56, -1, 137.3277     );
const ElemInfo LA("LA", 57, -1, 138.905477   );
const ElemInfo CE("CE", 58, -1, 140.1161     );
const ElemInfo PR("PR", 59, -1, 140.907652   );
const ElemInfo ND("ND", 60, -1, 144.2423     );
const ElemInfo PM("PM", 61, -1, 145.0        );
const ElemInfo SM("SM", 62, -1, 150.362      );
const ElemInfo EU("EU", 63, -1, 151.9641     );
const ElemInfo GD("GD", 64, -1, 157.253      );
const ElemInfo TB("TB", 65, -1, 158.925352   );
const ElemInfo DY("DY", 66, -1, 162.5001     );
const ElemInfo HO("HO", 67, -1, 164.930322   );
const ElemInfo ER("ER", 68, -1, 167.2593     );
const ElemInfo TM("TM", 69, -1, 168.934212   );
const ElemInfo YB("YB", 70, -1, 173.043      );
const ElemInfo LU("LU", 71, -1, 174.9671     );
const ElemInfo HF("HF", 72, -1, 178.492      );
const ElemInfo TA("TA", 73, -1, 180.947882   );
const ElemInfo W ("W" , 74, -1, 183.841      );
const ElemInfo RE("RE", 75, -1, 186.2071     );
const ElemInfo OS("OS", 76, -1, 190.233      );
const ElemInfo IR("IR", 77, -1, 192.2173     );
const ElemInfo PT("PT", 78, -1, 195.0849     );
const ElemInfo AU("AU", 79, -1, 196.9665694  );
const ElemInfo HG("HG", 80, -1, 200.592      );
const ElemInfo TL("TL", 81, -1, 204.38332    );
const ElemInfo PB("PB", 82, -1, 207.21       );
const ElemInfo BI("BI", 83, -1, 208.980401   );
const ElemInfo PO("PO", 84, -1, 210.0        );
const ElemInfo AT("AT", 85, -1, 210.0        );
const ElemInfo RN("RN", 86, -1, 220.0        );
const ElemInfo FR("FR", 87, -1, 223.0        );
const ElemInfo RA("RA", 88, -1, 226.0        );
const ElemInfo AC("AC", 89, -1, 227.0        );
const ElemInfo TH("TH", 90, -1, 232.038062   );
const ElemInfo PA("PA", 91, -1, 231.035882   );
const ElemInfo U ("U" , 92, -1, 238.028913   );
const ElemInfo NP("NP", 93, -1, 237.0        );
const ElemInfo PU("PU", 94, -1, 244.0        );
const ElemInfo AM("AM", 95, -1, 243.0        );
const ElemInfo CM("CM", 96, -1, 247.0        );
const ElemInfo BK("BK", 97, -1, 247.0        );
const ElemInfo CF("CF", 98, -1, 251.0        );
const ElemInfo ES("ES", 99, -1, 252.0        );
const ElemInfo FM("FM",100, -1, 257.0        );
const ElemInfo MD("MD",101, -1, 258.0        );
const ElemInfo NO("NO",102, -1, 259.0        );
const ElemInfo LR("LR",103, -1, 262.0        );
const ElemInfo RF("RF",104, -1, 261.0        );
const ElemInfo DB("DB",105, -1, 262.0        );
const ElemInfo SG("SG",106, -1, 266.0        );
const ElemInfo BH("BH",107, -1, 264.0        );
const ElemInfo HS("HS",108, -1, 277.0        );
const ElemInfo MT("MT",109, -1, 268.0        );
const ElemInfo DS("DS",110, -1, 271.0        );
const ElemInfo RG("RG",111, -1, 272.0        );
const ElemInfo CN("CN",112, -1, 285.0        );
const ElemInfo NH("NH",113, -1, 284.0        );
const ElemInfo FL("FL",114, -1, 289.0        );
const ElemInfo MC("MC",115, -1, 288.0        );
const ElemInfo LV("LV",116, -1, 292.0        );
const ElemInfo TS("TS",117, -1, 291.0        );
const ElemInfo OG("OG",118, -1, 294.0        );

constexpr int max_number   = 118;                              // 0-th element is NU
constexpr int max_period   = 8;                                // 0-th period is {NU}
constexpr int num_in_row[] = {0, 2, 10, 18, 36, 54, 86, 118};  // "NU" is not a pure element

const std::array<ElemInfo, max_number+1> ElementList = {
/* 0-th */  NU,
/* 1-th */  H ,                                                                 HE,
/* 2-th */  LI, BE,                                         B , C , N , O , F , NE,
/* 3-th */  NA, MG,                                         AL, SI, P , S , CL, AR,
/* 4-th */  K , CA, SC, TI, V , CR, MN, FE, CO, NI, CU, ZN, GA, GE, AS, SE, BR, KR,
/* 5-th */  RB, SR, Y , ZR, NB, MO, TC, RU, RH, PD, AG, CD, IN, SN, SB, TE, I , XE,
/* 6-th */  CS, BA,
/*Lanthanide*/  LA, CE, PR, ND, PM, SM, EU, GD, TB, DY, HO, ER, TM, YB, 
/* 6-th */          LU, HF, TA, W , RE, OS, IR, PT, AU, HG, TL, PB, BI, PO, AT, RN,
/* 7-th */  FR, RA,
/*Actinicles*/  AC, TH, PA, U , NP, PU, AM, CM, BK, CF, ES, FM, MD, NO,
/* 7-th */          LR, RF, DB, SG, BH, HS, MT, DS, RG, CN, NH, FL, MC, LV, TS, OG
};
// clang-format on

};  // namespace elem

static inline std::string getElemLabel(int Z) {
    if (Z <= 0 && Z > elem::max_number) throw std::runtime_error("bad Z number");
    return elem::ElementList[Z].label;
}

static inline int getElemIndex(const std::string& label) {
    std::string copylabel = label;
    std::for_each(copylabel.begin(), copylabel.end(), [](char& c) { c = ::toupper(c); });
    for (int z = 1; z <= elem::max_number; ++z) {
        if (copylabel == elem::ElementList[z].label) return z;
    }
    return 0;
}

static inline double getElemMass(int Z) {
    if (Z <= 0 && Z > elem::max_number) throw std::runtime_error("bad Z number");
    return elem::ElementList[Z].mass;
}
};  // namespace chem


// clang-format off
// const char* const ELEMENTS_LABEL[] = {
//  "NU" ,
//  "H" ,                                                                                "HE",
//  "LI","BE",                                                  "B" ,"C" ,"N" ,"O" ,"F" ,"NE",
//  "NA","MG",                                                  "AL","SI","P" ,"S" ,"CL","AR",
//  "K" ,"CA","SC","TI","V" ,"CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR",
//  "RB","SR","Y" ,"ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","I" ,"XE",
//  "CS","BA",
//       "LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB", // Lanthanide
//            "LU","HF","TA","W" ,"RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN",
//  "FR","RA",
//       "AC","TH","PA","U" ,"NP","PU","AM","CM","BK","CF","ES","FM","MD","NO", // Actinicles
//            "LR","RF","DB","SG","BH","HS","MT","DS","RG","CN","NH","FL","MC","LV","TS","Og"
// };

// const double ELEMENTS_MASS[] = // average of isotopic atoms
// {
// 0.0,            1.007947,       4.0026022,      6.9412,         9.0121823,
// 10.8117,        12.01078,       14.00672,       15.99943,       18.99840325,
// 20.17976,       22.989769282,   24.30506,       26.98153868,    28.08553,
// 30.9737622,     32.0655,        35.4532,        39.9481,        39.09831,
// 40.0784,        44.9559126,     47.8671,        50.94151,       51.99616,
// 54.9380455,     55.8452,        58.9331955,     58.69342,       63.5463,
// 65.4094,        69.7231,        72.641,         74.921602,      78.963,
// 79.9041,        83.7982,        85.46783,       87.621,         88.905852,
// 91.2242,        92.906382,      95.942,         98.0,           101.072,
// 102.905502,     106.421,        107.86822,      112.4118,       114.8183,
// 118.7107,       121.7601,       127.603,        126.904473,     131.2936,
// 132.90545192,   137.3277,       138.905477,     140.1161,       140.907652,
// 144.2423,       145.0,          150.362,        151.9641,       157.253,
// 158.925352,     162.5001,       164.930322,     167.2593,       168.934212,
// 173.043,        174.9671,       178.492,        180.947882,     183.841,
// 186.2071,       190.233,        192.2173,       195.0849,       196.9665694,
// 200.592,        204.38332,      207.21,         208.980401,     210.0,
// 210.0,          220.0,          223.0,          226.0,          227.0,
// 232.038062,     231.035882,     238.028913,     237.0,          244.0,
// 243.0,          247.0,          247.0,          251.0,          252.0,
// 257.0,          258.0,          259.0,          262.0,          261.0,
// 262.0,          266.0,          264.0,          277.0,          268.0,
// 271.0,          272.0,          285.0,          284.0,          289.0,
// 288.0,          292.0,          291.0,          294.0
// };

// const double ELEMENTS_MASS_NOAVG[] =
// {
// 0.0,                1.00782503207,      4.00260325415,      7.016004548,        9.012182201,
// 11.009305406,       12.0000000000,      14.00307400478,     15.99491461956,     18.998403224,
// 19.99244017542,     22.98976928087,     23.985041699,       26.981538627,       27.97692653246,
// 30.973761629,       31.972070999,       34.968852682,       39.96238312251,     38.963706679,
// 39.962590983,       44.955911909,       47.947946281,       50.943959507,       51.940507472,
// 54.938045141,       55.934937475,       58.933195048,       57.935342907,       62.929597474,
// 63.929142222,       68.925573587,       73.921177767,       74.921596478,       79.916521271,
// 78.918337087,       85.910610729,       84.911789737,       87.905612124,       88.905848295,
// 89.904704416,       92.906378058,       97.905408169,       98.906254747,       101.904349312,
// 102.905504292,      105.903485715,      106.90509682,       113.90335854,       114.903878484,
// 119.902194676,      120.903815686,      129.906224399,      126.904472681,      131.904153457,
// 132.905451932,      137.905247237,      138.906353267,      139.905438706,      140.907652769,
// 141.907723297,      144.912749023,      151.919732425,      152.921230339,      157.924103912,
// 158.925346757,      163.929174751,      164.93032207,       165.930293061,      168.93421325,
// 173.938862089,      174.940771819,      179.946549953,      180.947995763,      183.950931188,
// 186.955753109,      191.96148069,       192.96292643,       194.964791134,      196.966568662,
// 201.970643011,      204.974427541,      207.976652071,      208.980398734,      208.982430435,
// 210.987496271,      222.017577738,      222.01755173,       228.031070292,      227.027752127,
// 232.038055325,      231.03588399,       238.050788247,      237.048173444,      242.058742611,
// 243.06138108,       247.07035354,       247.07030708,       251.079586788,      252.082978512,
// 257.095104724,      258.098431319,      255.093241131,      260.105504,         263.112547,
// 255.107398,         259.114500,         262.122892,         263.128558,         265.136151,
// 281.162061,         272.153615,         283.171792,         283.176451,         285.183698,
// 287.191186,         292.199786,         291.206564,         293.214670
// };
// clang-format on

// namespace Elements {

// inline int ELEMENTS_ZNUM(const std::string& label) {
//     std::string copylabel = label;
//     std::for_each(copylabel.begin(), copylabel.end(), [](char& c) { c = ::toupper(c); });
//     for (int i = 1; i <= MAX_ELEMENTS_NUMBER; ++i) {
//         if (strcmp(copylabel.c_str(), ELEMENTS_LABEL[i]) == 0) return i;
//     }
//     return 0;
// }

// inline std::string ELEMENTS_NAME(const int& znumber) {
//     return (znumber > 0 && znumber <= MAX_ELEMENTS_NUMBER) ? std::string(ELEMENTS_LABEL[znumber]) : "NULL";
// }

// };  // namespace Elements


#endif  // CHEM_H
