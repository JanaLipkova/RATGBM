/*
 *  RATGBM.cp
 *  RATGBM
 *
 *  Created by Lipkova on 08/02/18.
 *  Copyright (c) 2018 Lipkova. All rights reserved.
 *
 */

#include "RATGBM.h"
#include "RATGBMPriv.h"

CFStringRef RATGBMUUID(void)
{
	CRATGBM* theObj = new CRATGBM;
	return theObj->UUID();
}

CFStringRef CRATGBM::UUID()
{
	return CFSTR("0001020304050607");
}
