//
// Julian.h
//
// Copyright (c) 2003-2012 Michael F. Henry
//
// Version 01/2013
//
#pragma once
#include "globals.h"
#include <string>

using namespace std;
using namespace OrbitTools;

// See note in Julian.cpp for information on this class and the epoch dates
const double EPOCH_JAN0_12H_1900 = 2415020.0; // Dec 31.5 1899 = Dec 31 1899 12h UTC
const double EPOCH_JAN1_00H_1900 = 2415020.5; // Jan  1.0 1900 = Jan  1 1900 00h UTC
const double EPOCH_JAN1_12H_2000 = 2451545.0; // Jan  1.5 2000 = Jan  1 2000 12h UTC

class Julian
{
public:
	Julian() { Initialize(2000, 1); }

	// Create a Julian date object from a year and day of year.
	// Example parameters: year = 2001, day = 1.5 (Jan 1 12h)
	Julian(int year, double day)
	{
		Initialize(year, day);
	}

	// Create a Julian date object.
	Julian(int year,               // i.e., 2004
		int mon,                // 1..12
		int day,                // 1..31
		int hour,               // 0..23
		int min,                // 0..59
		double sec /* = 0.0 */) // 0..(59.999999...)

	{
		// Calculate N, the day of the year (1..366)
		int N;
		int F1 = (int)((275.0 * mon) / 9.0);
		int F2 = (int)((mon + 9.0) / 12.0);

		if (IsLeapYear(year))
		{
			// Leap year
			N = F1 - F2 + day - 30;
		}
		else
		{
			// Common year
			N = F1 - (2 * F2) + day - 30;
		}

		double dblDay = N + (hour + (min + (sec / 60.0)) / 60.0) / 24.0;

		Initialize(year, dblDay);
	}

	~Julian() {};

	// Calculate Greenwich Mean Sidereal Time for the Julian date. The return value
	// is the angle, in radians, measuring eastward from the Vernal Equinox to the
	// prime meridian. This angle is also referred to as "ThetaG" (Theta GMST).
	// 
	// References:
	//    The 1992 Astronomical Almanac, page B6.
	//    Explanatory Supplement to the Astronomical Almanac, page 50.
	//    Orbital Coordinate Systems, Part III, Dr. T.S. Kelso, Satellite Times,
	//       Nov/Dec 1995
	double ToGmst() const
	{
		const double UT = fmod(m_Date + 0.5, 1.0);
		const double TU = (FromJan1_12h_2000() - UT) / 36525.0;

		double GMST = 24110.54841 + TU *
			(8640184.812866 + TU * (0.093104 - TU * 6.2e-06));

		GMST = fmod(GMST + SEC_PER_DAY * OMEGA_E * UT, SEC_PER_DAY);

		if (GMST < 0.0)
		{
			GMST += SEC_PER_DAY;  // "wrap" negative modulo value
		}

		return  (TWOPI * (GMST / SEC_PER_DAY));
	}

	// Calculate Local Mean Sidereal Time for given longitude (for this date).
	// The longitude is assumed to be in radians measured west from Greenwich.
	// The return value is the angle, in radians, measuring eastward from the
	// Vernal Equinox to the given longitude.
	double ToLmst(double lon) const
	{
		return fmod(ToGmst() + lon, TWOPI);
	}

	double FromJan0_12h_1900() const { return m_Date - EPOCH_JAN0_12H_1900; }
	double FromJan1_00h_1900() const { return m_Date - EPOCH_JAN1_00H_1900; }
	double FromJan1_12h_2000() const { return m_Date - EPOCH_JAN1_12H_2000; }

	// GetComponent()
	// Return requested components of date.
	// Year : Includes the century.
	// Month: 1..12
	// Day  : 1..31 including fractional part
	void GetComponent(int* pYear,
		int* pMon  /* = NULL */,
		double* pDOM  /* = NULL */) const
	{
		double jdAdj = Date() + 0.5;
		int    Z = (int)jdAdj;  // integer part
		double frac = jdAdj - Z;   // fractional part
		double alpha = trunc((Z - 1867216.25) / 36524.25);
		double A = Z + 1.0 + alpha - (int)(alpha / 4.0);
		double B = A + 1524.0;
		int    C = (int)((B - 122.1) / 365.25);
		int    D = (int)(C * 365.25);
		int    E = (int)((B - D) / 30.6001);

		double DOM = B - D - (int)(E * 30.6001) + frac;
		int    month = (E < 13.5) ? (E - 1) : (E - 13);
		int    year = (month > 2.5) ? (C - 4716) : (C - 4715);

		*pYear = year;

		if (pMon != NULL)
		{
			*pMon = month;
		}

		if (pDOM != NULL)
		{
			*pDOM = DOM;
		}
	}

	double Date() const { return m_Date; }

	void AddDay(double day) { m_Date += day; }
	void AddHour(double hr) { m_Date += (hr / HR_PER_DAY); }
	void AddMin(double min) { m_Date += (min / MIN_PER_DAY); }
	void AddSec(double sec) { m_Date += (sec / SEC_PER_DAY); }

	double SpanDay(const Julian& b) const { return m_Date - b.m_Date; }
	double SpanHour(const Julian& b) const { return SpanDay(b) * HR_PER_DAY; }
	double SpanMin(const Julian& b) const { return SpanDay(b) * MIN_PER_DAY; }
	double SpanSec(const Julian& b) const { return SpanDay(b) * SEC_PER_DAY; }

	static bool IsLeapYear(int y)
	{
		return (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
	}

protected:
	void Initialize(int year, double day)
	{
		// Calculate Julian date

		year--;

		// Centuries are not leap years unless they divide by 400
		int A = (year / 100);
		int B = 2 - A + (A / 4);

		double jan01 = trunc(365.25 * year) +
			trunc(30.6001 * 14) +
			1720994.5 + B;  // 1720994.5 = Oct 30, year -1

		m_Date = jan01 + day;
	}

	double m_Date; // Julian date
};