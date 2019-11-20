#include "pet_models.h"
#include "fwdmodel_linear1.h"

#include <fabber_core/fwdmodel.h>

extern "C" {
int CALL get_num_models()
{
    return 1;
}

const char *CALL get_model_name(int index)
{
    switch (index)
    {
    case 0:
        return "linear1";
        break;
    default:
        return NULL;
    }
}

NewInstanceFptr CALL get_new_instance_func(const char *name)
{
    if (string(name) == "linear1")
    {
        return PetFwdModel::NewInstance;
    }
    else
    {
        return NULL;
    }
}
}
