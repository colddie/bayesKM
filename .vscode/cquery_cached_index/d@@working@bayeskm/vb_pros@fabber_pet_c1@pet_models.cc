#include "pet_models.h"
#include "fwdmodel_pet_c1.h"

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
        return "pet";
        break;
    default:
        return NULL;
    }
}

NewInstanceFptr CALL get_new_instance_func(const char *name)
{
    if (string(name) == "pet")
    {
        return PetFwdModel::NewInstance;
    }
    else
    {
        return NULL;
    }
}
}
