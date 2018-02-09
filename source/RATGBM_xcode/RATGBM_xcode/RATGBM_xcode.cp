/*
 *  RATGBM_xcode.cp
 *  RATGBM_xcode
 *
 *  Created by Lipkova on 08/02/18.
 *  Copyright (c) 2018 Lipkova. All rights reserved.
 *
 */

#include "RATGBM_xcode.h"
#include "RATGBM_xcodePriv.h"

CFStringRef RATGBM_xcodeUUID(void)
{
	CRATGBM_xcode* theObj = new CRATGBM_xcode;
	return theObj->UUID();
}

CFStringRef CRATGBM_xcode::UUID()
{
	return CFSTR("0001020304050607");
}
