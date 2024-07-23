#ifndef MODEL_QMMMInterface_H
#define MODEL_QMMMInterface_H

class Model_QMMMInterface final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_QMMMInterface(){};

   private:
};

#endif  // MODEL_QMMMInterface_H