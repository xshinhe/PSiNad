#ifndef KIDS_PLATFORM_H
#define KIDS_PLATFORM_H

namespace PROJECT_NS {

class Platform : public std::enable_shared_from_this<Platform> {
   public:
    std::string getName();

   private:
    std::string _platname;
};


};  // namespace PROJECT_NS


#endif  // KIDS_PLATFORM_H
