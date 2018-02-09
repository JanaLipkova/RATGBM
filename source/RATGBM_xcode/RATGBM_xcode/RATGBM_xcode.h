/*
 *  RATGBM_xcode.h
 *  RATGBM_xcode
 *
 *  Created by Lipkova on 08/02/18.
 *  Copyright (c) 2018 Lipkova. All rights reserved.
 *
 */

extern "C" {
#include <CoreFoundation/CoreFoundation.h>

#pragma GCC visibility push(default)

/* External interface to the RATGBM_xcode, C-based */

CFStringRef RATGBM_xcodeUUID(void);

#pragma GCC visibility pop
}
