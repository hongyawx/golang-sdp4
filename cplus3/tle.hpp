#pragma once

#include <string>
#include "globals.h"
#include <map>

#include "math.h"
#include "Julian.hpp"
#include "Vector.hpp"

using namespace std;

class TLE
{
public:
	// Name
	const int TLE_LEN_LINE_DATA = 69; const int TLE_LEN_LINE_NAME = 24;

	// Line 1
	const int TLE1_COL_SATNUM = 2; const int TLE1_LEN_SATNUM = 5;
	const int TLE1_COL_INTLDESG_A = 9; const int TLE1_LEN_INTLDESG_A = 2;
	const int TLE1_COL_INTLDESG_B = 11; const int TLE1_LEN_INTLDESG_B = 3;
	const int TLE1_COL_INTLDESG_C = 14; const int TLE1_LEN_INTLDESG_C = 3;
	const int TLE1_COL_EPOCH_A = 18; const int TLE1_LEN_EPOCH_A = 2;
	const int TLE1_COL_EPOCH_B = 20; const int TLE1_LEN_EPOCH_B = 12;
	const int TLE1_COL_MEANMOTIONDT = 33; const int TLE1_LEN_MEANMOTIONDT = 10;
	const int TLE1_COL_MEANMOTIONDT2 = 44; const int TLE1_LEN_MEANMOTIONDT2 = 8;
	const int TLE1_COL_BSTAR = 53; const int TLE1_LEN_BSTAR = 8;
	const int TLE1_COL_EPHEMTYPE = 62; const int TLE1_LEN_EPHEMTYPE = 1;
	const int TLE1_COL_ELNUM = 64; const int TLE1_LEN_ELNUM = 4;

	// Line 2
	const int TLE2_COL_SATNUM = 2;  const int TLE2_LEN_SATNUM = 5;
	const int TLE2_COL_INCLINATION = 8;  const int TLE2_LEN_INCLINATION = 8;
	const int TLE2_COL_RAASCENDNODE = 17; const int TLE2_LEN_RAASCENDNODE = 8;
	const int TLE2_COL_ECCENTRICITY = 26; const int TLE2_LEN_ECCENTRICITY = 7;
	const int TLE2_COL_ARGPERIGEE = 34; const int TLE2_LEN_ARGPERIGEE = 8;
	const int TLE2_COL_MEANANOMALY = 43; const int TLE2_LEN_MEANANOMALY = 8;
	const int TLE2_COL_MEANMOTION = 52; const int TLE2_LEN_MEANMOTION = 11;
	const int TLE2_COL_REVATEPOCH = 63; const int TLE2_LEN_REVATEPOCH = 5;

	enum eTleLine
	{
		LINE_ZERO,
		LINE_ONE,
		LINE_TWO
	};

	enum eField
	{
		FLD_FIRST,
		FLD_NORADNUM = FLD_FIRST,
		FLD_INTLDESG,  // International designator
		FLD_SET,       // TLE set number
		FLD_EPOCHYEAR, // Epoch: Last two digits of year
		FLD_EPOCHDAY,  // Epoch: Fractional Julian Day of year
		FLD_ORBITNUM,  // Orbit at epoch
		FLD_I,         // Inclination
		FLD_RAAN,      // R.A. ascending node
		FLD_E,         // Eccentricity
		FLD_ARGPER,    // Argument of perigee
		FLD_M,         // Mean anomaly
		FLD_MMOTION,   // Mean motion
		FLD_MMOTIONDT, // First time derivative of mean motion
		FLD_MMOTIONDT2,// Second time derivative of mean motion
		FLD_BSTAR,     // BSTAR Drag
		FLD_LAST       // MUST be last
	};

	enum eUnits
	{
		U_FIRST,
		U_RAD = U_FIRST,  // radians
		U_DEG,            // degrees
		U_NATIVE,         // TLE format native units (no conversion)
		U_LAST            // MUST be last
	};
	
public:
	TLE(const string& strName, const string& strLine1, const string& strLine2)
	{
		m_strLine0 = strName;
		m_strLine1 = strLine1;
		m_strLine2 = strLine2;

		TrimRight(m_strLine0);
		Initialize();
	}

	bool model_init()
	{
		bool sgp_mode = false;
		InitializeCachingVars();

		int    epochYear = (int)GetField(FLD_EPOCHYEAR);
		double epochDay = GetField(FLD_EPOCHDAY);

		if (epochYear < 57)
		{
			epochYear += 2000;
		}
		else
		{
			epochYear += 1900;
		}

		m_jdEpoch = Julian(epochYear, epochDay);

		m_secPeriod = -1.0;

		// Recover the original mean motion and semimajor axis from the
		// input elements.
		double mm = MeanMotionTle();
		double rpmin = mm * TWOPI / MIN_PER_DAY;   // rads per minute

		double a1 = pow(XKE / rpmin, 2.0 / 3.0);
		double e = Eccentricity();
		double i = m_Inclination;
		double temp = (1.5 * CK2 * (3.0 * sqr(cos(i)) - 1.0) /
			pow(1.0 - e * e, 1.5));
		double delta1 = temp / (a1 * a1);
		double a0 = a1 *
			(1.0 - delta1 *
				((1.0 / 3.0) + delta1 *
					(1.0 + 134.0 / 81.0 * delta1)));

		double delta0 = temp / (a0 * a0);

		m_rmMeanMotionRec = rpmin / (1.0 + delta0);
		m_aeAxisSemiMajorRec = a0 / (1.0 - delta0);
		m_aeAxisSemiMinorRec = m_aeAxisSemiMajorRec * sqrt(1.0 - (e * e));
		m_kmPerigeeRec = XKMPER_WGS72 * (m_aeAxisSemiMajorRec * (1.0 - e) - AE);
		m_kmApogeeRec = XKMPER_WGS72 * (m_aeAxisSemiMajorRec * (1.0 + e) - AE);

		if (TWOPI / m_rmMeanMotionRec >= 225.0)
		{
			// SDP4 - period >= 225 minutes.
			//m_pNoradModel = new cNoradSDP4(*this);
			sgp_mode = false;
		}
		else
		{
			// SGP4 - period < 225 minutes
			//m_pNoradModel = new cNoradSGP4(*this);
			sgp_mode = true;
		}

		return sgp_mode;
	}

public:
	void Initialize()
	{
		// Have we already been initialized?
		if (m_Field[FLD_NORADNUM].size()) { return; }

		m_Field[FLD_NORADNUM] = m_strLine1.substr(TLE1_COL_SATNUM, TLE1_LEN_SATNUM);
		m_Field[FLD_INTLDESG] = m_strLine1.substr(TLE1_COL_INTLDESG_A,
			TLE1_LEN_INTLDESG_A +
			TLE1_LEN_INTLDESG_B +
			TLE1_LEN_INTLDESG_C);
		m_Field[FLD_EPOCHYEAR] =
			m_strLine1.substr(TLE1_COL_EPOCH_A, TLE1_LEN_EPOCH_A);

		m_Field[FLD_EPOCHDAY] =
			m_strLine1.substr(TLE1_COL_EPOCH_B, TLE1_LEN_EPOCH_B);

		if (m_strLine1[TLE1_COL_MEANMOTIONDT] == '-')
		{
			// value is negative
			m_Field[FLD_MMOTIONDT] = "-0";
		}
		else
		{
			m_Field[FLD_MMOTIONDT] = "0";
		}

		m_Field[FLD_MMOTIONDT] += m_strLine1.substr(TLE1_COL_MEANMOTIONDT + 1,
			TLE1_LEN_MEANMOTIONDT);

		// decimal point assumed; exponential notation
		m_Field[FLD_MMOTIONDT2] = ExpToAtof(
			m_strLine1.substr(TLE1_COL_MEANMOTIONDT2,
				TLE1_LEN_MEANMOTIONDT2));
		// decimal point assumed; exponential notation
		m_Field[FLD_BSTAR] = ExpToAtof(m_strLine1.substr(TLE1_COL_BSTAR,
			TLE1_LEN_BSTAR));
		// TLE1_COL_EPHEMTYPE      
		// TLE1_LEN_EPHEMTYPE   
		m_Field[FLD_SET] = m_strLine1.substr(TLE1_COL_ELNUM, TLE1_LEN_ELNUM);

		TrimLeft(m_Field[FLD_SET]);

		// TLE2_COL_SATNUM         
		// TLE2_LEN_SATNUM         

		m_Field[FLD_I] = m_strLine2.substr(TLE2_COL_INCLINATION,
			TLE2_LEN_INCLINATION);
		TrimLeft(m_Field[FLD_I]);

		m_Field[FLD_RAAN] = m_strLine2.substr(TLE2_COL_RAASCENDNODE,
			TLE2_LEN_RAASCENDNODE);
		TrimLeft(m_Field[FLD_RAAN]);

		// decimal point is assumed
		m_Field[FLD_E] = "0.";
		m_Field[FLD_E] += m_strLine2.substr(TLE2_COL_ECCENTRICITY,
			TLE2_LEN_ECCENTRICITY);

		m_Field[FLD_ARGPER] = m_strLine2.substr(TLE2_COL_ARGPERIGEE,
			TLE2_LEN_ARGPERIGEE);
		TrimLeft(m_Field[FLD_ARGPER]);

		m_Field[FLD_M] = m_strLine2.substr(TLE2_COL_MEANANOMALY,
			TLE2_LEN_MEANANOMALY);
		TrimLeft(m_Field[FLD_M]);

		m_Field[FLD_MMOTION] = m_strLine2.substr(TLE2_COL_MEANMOTION,
			TLE2_LEN_MEANMOTION);
		TrimLeft(m_Field[FLD_MMOTION]);

		m_Field[FLD_ORBITNUM] = m_strLine2.substr(TLE2_COL_REVATEPOCH,
			TLE2_LEN_REVATEPOCH);
		TrimLeft(m_Field[FLD_ORBITNUM]);
	}

	static string GetUnits(eField fld)
	{
		static const string strDegrees = " degrees";
		static const string strRevsPerDay = " revs / day";
		static const string strNull;

		switch (fld)
		{
		case FLD_I:
		case FLD_RAAN:
		case FLD_ARGPER:
		case FLD_M:
			return strDegrees;

		case FLD_MMOTION:
			return strRevsPerDay;

		default:
			return strNull;
		}
	}

	static double ConvertUnits(double valNative, eField fld, eUnits units)
	{
		switch (fld)
		{
		case FLD_I:
		case FLD_RAAN:
		case FLD_ARGPER:
		case FLD_M:
		{
			// The native TLE format is DEGREES
			if (units == U_RAD)
			{
				return valNative * RADS_PER_DEG;
			}
		}
		}

		return valNative; // return value in unconverted native format
	}

	//  Return requested field as a double (function return value) or as a text
	double GetField(eField fld,  // which field to retrieve
		eUnits units = U_NATIVE,  // return units in rad, deg etc.
		string* pstr = NULL,      // return ptr for str value
		bool bStrUnits = false)   // 'true': append units to str val
		const
	{
		if (pstr)
		{
			// Return requested field in string form.
			*pstr = m_Field[fld];

			if (bStrUnits)
			{
				*pstr += GetUnits(fld);
			}

			TrimLeft(*pstr);
			TrimRight(*pstr);

			return 0.0;
		}
		else
		{
			// Return requested field in floating-point form.
			// Return cache contents if it exists, else populate cache
			FldKey key = Key(units, fld);

			if (m_mapCache.find(key) == m_mapCache.end())
			{
				// Value not in cache; add it
				double valNative = atof(m_Field[fld].c_str());
				double valConv = ConvertUnits(valNative, fld, units);
				m_mapCache[key] = valConv;

				return valConv;
			}
			else
			{
				// return cached value
				return m_mapCache[key];
			}
		}
	}

	static string ExpToAtof(const string& str)
	{
		const int COL_SIGN = 0;
		const int LEN_SIGN = 1;

		const int COL_MANTISSA = 1;
		const int LEN_MANTISSA = 5;

		const int COL_EXPONENT = 6;
		const int LEN_EXPONENT = 2;

		string sign = str.substr(COL_SIGN, LEN_SIGN);
		string mantissa = str.substr(COL_MANTISSA, LEN_MANTISSA);
		string exponent = str.substr(COL_EXPONENT, LEN_EXPONENT);

		TrimLeft(exponent);

		return sign + "0." + mantissa + "e" + exponent;
	}

	static void TrimLeft(string& s)
	{
		while ((s.size() > 0) && (s[0] == ' '))
		{
			s.erase(0, 1);
		}
	}

	static void TrimRight(string& s)
	{
		while ((s.size() > 0) && (s[s.size() - 1] == ' '))
		{
			s.erase(s.size() - 1);
		}
	}

	double RadGet(eField fld) const { return GetField(fld, U_RAD); }
	double DegGet(eField fld) const { return GetField(fld, U_DEG); }

	double Inclination()   const { return m_Inclination; }
	double Eccentricity()  const { return m_Eccentricity; }
	double RAAN()          const { return m_RAAN; }
	double ArgPerigee()    const { return m_ArgPerigee; }
	double BStar()         const { return m_BStar; }
	double MeanMotionTle() const { return m_TleMeanMotion; }
	double MeanAnomaly()   const { return m_MeanAnomaly; }

	// "Recovered" from the input elements
	double SemiMajor()  const { return m_aeAxisSemiMajorRec; }
	double MeanMotion() const { return m_rmMeanMotionRec; }  // mean motion, rads/min

	Julian Epoch() const { return m_jdEpoch; }

	void InitializeCachingVars()
	{
		m_Inclination = RadGet(FLD_I);
		m_Eccentricity = GetField(FLD_E);
		m_RAAN = RadGet(FLD_RAAN);
		m_ArgPerigee = RadGet(FLD_ARGPER);
		m_BStar = GetField(FLD_BSTAR) / AE;
		//m_Drag = GetField(FLD_MMOTIONDT);
		m_TleMeanMotion = GetField(FLD_MMOTION);
		m_MeanAnomaly = RadGet(FLD_M);
		//m_RevAtEpoch = (int)GetField(FLD_ORBITNUM);
		//m_TleSetNumber = (int)GetField(FLD_SET);
	}

public:
	// Satellite name and two data lines
	string m_strLine0;
	string m_strLine1;
	string m_strLine2;

	// Converted fields, in atof()-readable form
	string m_Field[FLD_LAST];

	// Cache of field values in "double" format
	typedef int FldKey;
	FldKey Key(eUnits u, eField f) const { return (u * 100) + f; }
	mutable map<FldKey, double>  m_mapCache;

	// Caching variables; note units are not necessarily the same as TLE units
	mutable double m_secPeriod;

	// Caching variables for standard TLE elements
	double m_Inclination;
	double m_Eccentricity;
	double m_RAAN;
	double m_ArgPerigee;
	double m_BStar;
	double m_TleMeanMotion;
	double m_MeanAnomaly;

	// Caching variables recovered from the input TLE elements
	double m_aeAxisSemiMajorRec;  // semi-major axis, in AE units
	double m_aeAxisSemiMinorRec;  // semi-minor axis, in AE units
	double m_rmMeanMotionRec;     // radians per minute
	double m_kmPerigeeRec;        // perigee, in km
	double m_kmApogeeRec;         // apogee, in km

	// Julian date
	Julian m_jdEpoch;
};