#ifndef __PANDA_SUBSYSTEM_H__
#define __PANDA_SUBSYSTEM_H__

#include "TObject.h"
#include "string.h"
using namespace std;
class PandaSubsystem : public TObject {

    public:

    std::string subsystem_name;

    PandaSubsystem();
    ~PandaSubsystem();

	ClassDef(PandaSubsystem, 1)
};

#endif
