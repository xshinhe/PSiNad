#ifndef OFSManager_H
#define OFSManager_H

#include <fstream>
#include <iostream>
#include <map>
#include <vector>

class OFSManager {
   public:
    OFSManager(){};

    int push(const std::string& str) {
        if (ofs_ptr_loc.find(str) != ofs_ptr_loc.end()) return ofs_ptr_loc[str];
        ofs_ptr_loc[str] = ofs_ptr_arr.size();
        ofs_ptr_arr.push_back(new std::ofstream(str));
        return ofs_ptr_loc[str];
    }

    inline std::ofstream& stream(int id) { return *ofs_ptr_arr[id]; }

    inline std::ofstream& stream(std::string key) { return *ofs_ptr_arr[ofs_ptr_loc[key]]; }

    virtual ~OFSManager() {
        for (auto&& ofs_ptr : ofs_ptr_arr) {
            ofs_ptr->close();
            delete ofs_ptr;
        }
    }

    std::map<std::string, int> ofs_ptr_loc;
    std::vector<std::ofstream*> ofs_ptr_arr;
};


#endif  // OFSManager_H
