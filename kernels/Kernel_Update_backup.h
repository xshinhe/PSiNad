#ifndef Kernel_Update_backup_H
#define Kernel_Update_backup_H

class Kernel_Update_backup final : public Kernel {
   private:
    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

#endif  // Kernel_Update_backup_H
