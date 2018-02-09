/*
 *  RATGBM.h
 *  RATGBM
 *
 *  Created by Lipkova on 08/02/18.
 *  Copyright (c) 2018 Lipkova. All rights reserved.
 *
 */

extern "C" {
#include <CoreFoundation/CoreFoundation.h>

#pragma GCC visibility push(default)

/* External interface to the RATGBM, C-based */

CFStringRef RATGBMUUID(void);

#pragma GCC visibility pop
}
