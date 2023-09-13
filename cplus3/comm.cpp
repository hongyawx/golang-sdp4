#include <iostream>
#include <cstring>
#include "geo.hpp"
#include "SGP4.hpp"
#include "SDP4.hpp"
#include "comm.h"

using namespace std;

struct AngleResult calculateAngle(const char*  line1, const char*  line2, const char*  line3, int year, int month, int day, int hour, int min, int sec, double lat, double lon, double altKm)
{
	struct AngleResult result = {};

    std::string line11(line1);
    std::string line22(line2);
    std::string line33(line3);

	TLE tle(line11, line22, line33);
	bool sgp_model = tle.model_init();
	if (sgp_model)
	{
		cout << "初始化成功" << endl;
	}
	else
	{
		cout << "初始化失败" << endl;
	}
	NoradBase* norad_base = new SGP4(tle);
	Julian ju_now(year, month, day, hour, min, sec);
	double mpe = ju_now.SpanMin(tle.Epoch());
	norad_base->GetPosition(mpe);
	double radiusAe = XKMPER_WGS72 / AE;
	norad_base->ScalePosVector(radiusAe);
	norad_base->ScaleVelVector(radiusAe * (MIN_PER_DAY / 86400));
	Geo geo(lat, lon, altKm);
	geo.GetLookAngle(norad_base->Position(), norad_base->Velocity(), norad_base->Date());
	result.pitch = geo.ElevationDeg();
	result.yaw = geo.AzimuthDeg();
	delete norad_base;
	return result;
}
