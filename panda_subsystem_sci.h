#ifndef __PANDA_SUBSYSTEM_SCI_H__
#define __PANDA_SUBSYSTEM_SCI_H__

#include "SciEvent.h"
#include "panda_subsystem.h"
#include "string.h"

class PandaSubsystemSCI : public PandaSubsystem
{

  public:
    SciEvent sci_raw;

    PandaSubsystemSCI();
    ~PandaSubsystemSCI();

    ClassDef(PandaSubsystemSCI, 1)
};

#endif
