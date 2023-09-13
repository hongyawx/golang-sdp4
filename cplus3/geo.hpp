#pragma once

#include <string>
#include "globals.h"

#include "Julian.hpp"
#include "Vector.hpp"

// Geocentric coordinates.
class Geo
{
public:
	Geo(double latRad, double lonRad, double altKm)
		: m_Lat(deg2rad(latRad)),
		m_Lon(deg2rad(lonRad)),
		m_Alt(altKm)
	{
	}

	virtual ~Geo() {}

	double LatitudeRad()  const { return m_Lat; }
	double LongitudeRad() const { return m_Lon; }

	double LatitudeDeg()  const { return rad2deg(m_Lat); }
	double LongitudeDeg() const { return rad2deg(m_Lon); }

	double AltitudeKm() const { return m_Alt; }
	void   AltitudeKm(double alt) { m_Alt = alt; }

	const Vector& Position() const { return m_Position; }
	const Vector& Velocity() const { return m_Velocity; }

	double AzimuthRad()     const { return m_Az; }
	double ElevationRad()   const { return m_El; }

	double AzimuthDeg()     const { return rad2deg(m_Az); }
	double ElevationDeg()   const { return rad2deg(m_El); }

	double RangeKm()        const { return m_Range; }
	double RangeRateKmSec() const { return m_RangeRate; }

	// Creates a instance of the class from geodetic coordinates.
	//
	// Assumes the Earth is an oblate spheroid.
	// Reference: The 1992 Astronomical Almanac, page K11
	// Reference: www.celestrak.com (Dr. T.S. Kelso)
	void eci_posi(Julian date)
	{
		double lat = LatitudeRad();
		double lon = LongitudeRad();
		double alt = AltitudeKm();

		// Calculate Local Mean Sidereal Time (theta)
		double theta = date.ToLmst(lon);
		double c = 1.0 / sqrt(1.0 + F * (F - 2.0) * sqr(sin(lat)));
		double s = sqr(1.0 - F) * c;
		double achcp = (XKMPER_WGS72 * c + alt) * cos(lat);

		m_Position.m_x = achcp * cos(theta);         // km
		m_Position.m_y = achcp * sin(theta);         // km
		m_Position.m_z = (XKMPER_WGS72 * s + alt) * sin(lat);   // km
		m_Position.m_w = sqrt(sqr(m_Position.m_x) +  // range, km
			sqr(m_Position.m_y) +
			sqr(m_Position.m_z));

		// Determine velocity components due to earth's rotation
		double mfactor = TWOPI * (OMEGA_E / SEC_PER_DAY);

		m_Velocity.m_x = -mfactor * m_Position.m_y;  // km / sec
		m_Velocity.m_y = mfactor * m_Position.m_x;  // km / sec
		m_Velocity.m_z = 0.0;                        // km / sec
		m_Velocity.m_w = sqrt(sqr(m_Velocity.m_x) +  // range rate km/sec^2
			sqr(m_Velocity.m_y));
	}

	// Return the topocentric (azimuth, elevation, etc.) coordinates for a target
	// object located at the given ECI coordinates.
	void GetLookAngle(Vector pos, Vector vel, Julian date)
	{
		// Calculate the ECI coordinates for this cSite object at the time of interest.
		eci_posi(date);

		Vector vecRgRate(vel.m_x - Velocity().m_x,
			vel.m_y - Velocity().m_y,
			vel.m_z - Velocity().m_z);

		double x = pos.m_x - Position().m_x;
		double y = pos.m_y - Position().m_y;
		double z = pos.m_z - Position().m_z;
		double w = sqrt(sqr(x) + sqr(y) + sqr(z));

		Vector vecRange(x, y, z, w);

		// The site's Local Mean Sidereal Time at the time of interest.
		double theta = date.ToLmst(LongitudeRad());

		double sin_lat = sin(LatitudeRad());
		double cos_lat = cos(LatitudeRad());
		double sin_theta = sin(theta);
		double cos_theta = cos(theta);

		double top_s = sin_lat * cos_theta * vecRange.m_x +
			sin_lat * sin_theta * vecRange.m_y -
			cos_lat * vecRange.m_z;
		double top_e = -sin_theta * vecRange.m_x +
			cos_theta * vecRange.m_y;
		double top_z = cos_lat * cos_theta * vecRange.m_x +
			cos_lat * sin_theta * vecRange.m_y +
			sin_lat * vecRange.m_z;
		double az = atan(-top_e / top_s);

		if (top_s > 0.0)
		{
			az += PI;
		}

		if (az < 0.0)
		{
			az += 2.0 * PI;
		}

		double el = asin(top_z / vecRange.m_w);
		double rate = (vecRange.m_x * vecRgRate.m_x +
			vecRange.m_y * vecRgRate.m_y +
			vecRange.m_z * vecRgRate.m_z) / vecRange.m_w;

#ifdef WANT_ATMOSPHERIC_CORRECTION
		double saveEl = el;

		// Elevation correction for atmospheric refraction.
		// Reference:  Astronomical Algorithms by Jean Meeus, pp. 101-104
		// Note:  Correction is meaningless when apparent elevation is below horizon
		el += deg2rad((1.02 /
			tan(deg2rad(rad2deg(el) + 10.3 /
				(rad2deg(el) + 5.11)))) / 60.0);
		if (el < 0.0)
		{
			// Reset to true elevation
			el = saveEl;
		}

		if (el > (PI / 2))
		{
			el = (PI / 2);
		}
#endif

		m_Az = az;
		m_El = el;
		m_Range = vecRange.m_w;
		m_RangeRate = rate;
	}

public:
	double m_Lat;   // Latitude,  radians (negative south)
	double m_Lon;   // Longitude, radians (negative west)
	double m_Alt;   // Altitude,  km      (above ellipsoid height)

	Vector  m_Position;
	Vector  m_Velocity;

	double m_Az;         // Azimuth, radians
	double m_El;         // Elevation, radians
	double m_Range;      // Range, kilometers
	double m_RangeRate;  // Range rate of change, km/sec
						 // Negative value means "towards observer"
};