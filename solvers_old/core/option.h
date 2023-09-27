#ifndef Option_H
#define Option_H

// template<typename T>

struct Option {
    std::string flag;
    int type;
    std::map<std::string, int> _dict;
}

template <typename Op>
class OptionImpl {};



#endif  // Option_H
