#ifndef MODEL_QMInterface_H
#define MODEL_QMInterface_H

DEFINE_POLICY(QMPolicy,  //
              GAUSSIAN,  //
              ORCA,      //
              MNDO,      //
              BAGEL,     //
              MOLCAS,    //
              NONE);     //

class Model_QMInterface final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_QMInterface(){};

   private:
};

#endif  // MODEL_QMInterface_H