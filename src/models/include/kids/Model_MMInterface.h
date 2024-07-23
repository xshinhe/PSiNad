#ifndef MODEL_MMInterface_H
#define MODEL_MMInterface_H

#include "kids/Model.h"
#include "kids/Model_Bath.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(MMPolicy,  //
              Openmm,    //
              Amber,     //
              Gromacs,   //
              None);     //

class Model_MMInterface final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_MMInterface(){};

   private:
};

};  // namespace PROJECT_NS


#endif  // MODEL_MMInterface_H