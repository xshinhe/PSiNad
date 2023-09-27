/**
 *
 * This file provide basic abstract of class & interface
 *
 * All class in this program should be the child of this class.
 *
 */

#ifndef MODEL_H
#define MODEL_H

#include <string>

#include "../utils/param_wrapper.h"

struct Option {  // @TODO
    int type;
    std::string flag;
};

struct IOUnit {  // @TODO simplify code
    double leng = 1, time = 1, mass = 1, ener = 1, temp = 1;
};

/**
 * @brief Model is a basic abstract class and interface in simulations.
 * @detail in virtue of that everything (theorem & object) is no more than a
 * model. It provides nothing, other than storing a parameter list and set
 * universal units.
 */
class Model {
   public:
    /**
     * @brief constructor from Param
     */
    Model(const Param& iparm);

    /**
     * @brief constructor from string (Json-string). API for python
     */
    Model(const std::string& iparm_str) : Model(Param::parse(iparm_str)){};

    /**
     * @brief destructor (pure virtual make it an abstract interface)
     */
    virtual ~Model() = 0;

    /**
     * @brief static interface parsing unit parm to data structure
     */
    static IOUnit parse_unit_parm(const Param& parm);

    IOUnit iou;  ///< Units used for runtime IO

    std::string tag;  ///< public tag

   protected:
    DEFINE_POINTER_PROTECTED(num_real, workr);     ///< allocate real workspace if it's necessary
    DEFINE_POINTER_PROTECTED(num_complex, workc);  ///< allocate complex workspace if it's necessary

    Param parm;  ///< Param structure
};

#endif  // MODEL_H
