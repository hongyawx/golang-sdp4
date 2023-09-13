// comm.h
#ifdef __cplusplus
extern "C" {
#endif


struct AngleResult
{
	float pitch;
	float yaw;
};

struct AngleResult calculateAngle(const char*  line1, const char*  line2, const char*  line3, int year, int month, int day, int hour, int min, int sec, double lat, double lon, double altKm) ;


#ifdef __cplusplus
}
#endif