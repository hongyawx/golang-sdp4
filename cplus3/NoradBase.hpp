#pragma once

#include <string>
#include "tle.hpp"

class NoradBase
{
public:
	NoradBase(const TLE& orbit) :
		m_Orbit(orbit)
	{
		// Initialize any variables which are time-independent when
		// calculating the ECI coordinates of the satellite.
		m_sinio = sin(m_Orbit.Inclination());
		m_cosio = cos(m_Orbit.Inclination());

		double theta2 = m_cosio * m_cosio;
		double x3thm1 = 3.0 * theta2 - 1.0;
		double eosq = sqr(m_Orbit.Eccentricity());

		m_betao2 = 1.0 - eosq;
		m_betao = sqrt(m_betao2);

		// For perigee below 156 km, the values of S and QOMS2T are altered.
		double rp = m_Orbit.SemiMajor() * (1.0 - m_Orbit.Eccentricity());
		double perigee = (rp - 1.0) * XKMPER_WGS72;

		m_s4 = S;
		m_qoms24 = QOMS2T;

		if (perigee < 156.0)
		{
			m_s4 = perigee - 78.0;

			if (perigee <= 98.0)
			{
				m_s4 = 20.0;
			}

			m_qoms24 = pow((120.0 - m_s4) * AE / XKMPER_WGS72, 4.0);
			m_s4 = m_s4 / XKMPER_WGS72 + AE;
		}

		const double pinvsq = 1.0 / (sqr(m_Orbit.SemiMajor()) * sqr(m_betao2));

		m_tsi = 1.0 / (m_Orbit.SemiMajor() - m_s4);
		m_eta = m_Orbit.SemiMajor() * m_Orbit.Eccentricity() * m_tsi;
		m_eeta = m_Orbit.Eccentricity() * m_eta;

		const double etasq = m_eta * m_eta;
		const double psisq = fabs(1.0 - etasq);

		m_coef = m_qoms24 * pow(m_tsi, 4.0);
		m_coef1 = m_coef / pow(psisq, 3.5);

		const double c2 = m_coef1 * m_Orbit.MeanMotion() *
			(m_Orbit.SemiMajor() * (1.0 + 1.5 * etasq + m_eeta * (4.0 + etasq)) +
				0.75 * CK2 * m_tsi / psisq * x3thm1 *
				(8.0 + 3.0 * etasq * (8.0 + etasq)));

		m_c1 = m_Orbit.BStar() * c2;
		m_a3ovk2 = -XJ3 / CK2 * pow(AE, 3.0);

		m_c3 = m_coef * m_tsi * m_a3ovk2 * m_Orbit.MeanMotion() * AE * m_sinio / m_Orbit.Eccentricity();

		const double x1mth2 = 1.0 - theta2;
		m_c4 = 2.0 * m_Orbit.MeanMotion() * m_coef1 * m_Orbit.SemiMajor() * m_betao2 *
			(m_eta * (2.0 + 0.5 * etasq) +
				m_Orbit.Eccentricity() * (0.5 + 2.0 * etasq) -
				2.0 * CK2 * m_tsi / (m_Orbit.SemiMajor() * psisq) *
				(-3.0 * x3thm1 * (1.0 - 2.0 * m_eeta + etasq * (1.5 - 0.5 * m_eeta)) +
					0.75 * x1mth2 *
					(2.0 * etasq - m_eeta * (1.0 + etasq)) *
					cos(2.0 * m_Orbit.ArgPerigee())));

		const double theta4 = theta2 * theta2;
		const double temp1 = 3.0 * CK2 * pinvsq * m_Orbit.MeanMotion();;
		const double temp2 = temp1 * CK2 * pinvsq;
		const double temp3 = 1.25 * CK4 * pinvsq * pinvsq * m_Orbit.MeanMotion();;

		m_xmdot = m_Orbit.MeanMotion() + 0.5 * temp1 * m_betao * x3thm1 +
			0.0625 * temp2 * m_betao *
			(13.0 - 78.0 * theta2 + 137.0 * theta4);

		const double x1m5th = 1.0 - 5.0 * theta2;

		m_omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 *
			(7.0 - 114.0 * theta2 + 395.0 * theta4) +
			temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);

		const double xhdot1 = -temp1 * m_cosio;

		m_xnodot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) +
			2.0 * temp3 * (3.0 - 7.0 * theta2)) * m_cosio;
		m_xnodcf = 3.5 * m_betao2 * xhdot1 * m_c1;
		m_t2cof = 1.5 * m_c1;
	}

	virtual ~NoradBase() { }

	virtual void GetPosition(double tsince) = 0;

protected:
	void FinalPosition(double incl, double omega, double e, double a,
		double xl, double xnode, double xn, double tsince)
	{
		if ((e * e) > 1.0)
		{
			// Error in satellite data
			e = 0.1;
		}

		double beta = sqrt(1.0 - e * e);

		// Long period periodics 
		double axn = e * cos(omega);
		double temp = 1.0 / (a * beta * beta);

		double sinip = sin(m_Orbit.Inclination());
		double cosip = cos(m_Orbit.Inclination());
		double aycof = 0.25 * m_a3ovk2 * sinip;
		double xlcof = (0.125 * m_a3ovk2 * sinip * (3.0 + 5.0 * cosip)) /
			(1.0 + cosip);
		double xll = temp * xlcof * axn;
		double aynl = temp * aycof;
		double xlt = xl + xll;
		double ayn = e * sin(omega) + aynl;

		const double E6A = 1.0e-06;

		// Solve Kepler's Equation 

		double capu = Fmod2p(xlt - xnode);
		double temp2 = capu;
		double temp3 = 0.0;
		double temp4 = 0.0;
		double temp5 = 0.0;
		double temp6 = 0.0;
		double sinepw = 0.0;
		double cosepw = 0.0;
		bool   fDone = false;

		for (int i = 1; (i <= 10) && !fDone; i++)
		{
			sinepw = sin(temp2);
			cosepw = cos(temp2);
			temp3 = axn * sinepw;
			temp4 = ayn * cosepw;
			temp5 = axn * cosepw;
			temp6 = ayn * sinepw;

			double epw = (capu - temp4 + temp3 - temp2) /
				(1.0 - temp5 - temp6) + temp2;

			if (fabs(epw - temp2) <= E6A)
			{
				fDone = true;
			}
			else
			{
				temp2 = epw;
			}
		}

		// Short period preliminary quantities 
		double ecose = temp5 + temp6;
		double esine = temp3 - temp4;
		double elsq = axn * axn + ayn * ayn;
		temp = 1.0 - elsq;
		double pl = a * temp;
		double r = a * (1.0 - ecose);
		double temp1 = 1.0 / r;
		double rdot = XKE * sqrt(a) * esine * temp1;
		double rfdot = XKE * sqrt(pl) * temp1;
		temp2 = a * temp1;
		double betal = sqrt(temp);
		temp3 = 1.0 / (1.0 + betal);
		double cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
		double sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
		double u = AcTan(sinu, cosu);
		double sin2u = 2.0 * sinu * cosu;
		double cos2u = 2.0 * cosu * cosu - 1.0;

		temp = 1.0 / pl;
		temp1 = CK2 * temp;
		temp2 = temp1 * temp;

		// Update for short periodics 
		double cosip2 = cosip * cosip;
		double x3thm1 = 3.0 * cosip2 - 1.0;
		double x1mth2 = 1.0 - cosip2;
		double x7thm1 = 7.0 * cosip2 - 1.0;
		double rk = r * (1.0 - 1.5 * temp2 * betal * x3thm1) +
			0.5 * temp1 * x1mth2 * cos2u;
		double uk = u - 0.25 * temp2 * x7thm1 * sin2u;
		double xnodek = xnode + 1.5 * temp2 * m_cosio * sin2u;
		double xinck = incl + 1.5 * temp2 * m_cosio * m_sinio * cos2u;
		double rdotk = rdot - xn * temp1 * x1mth2 * sin2u;
		double rfdotk = rfdot + xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1);

		// Orientation vectors 
		double sinuk = sin(uk);
		double cosuk = cos(uk);
		double sinik = sin(xinck);
		double cosik = cos(xinck);
		double sinnok = sin(xnodek);
		double cosnok = cos(xnodek);
		double xmx = -sinnok * cosik;
		double xmy = cosnok * cosik;
		double ux = xmx * sinuk + cosnok * cosuk;
		double uy = xmy * sinuk + sinnok * cosuk;
		double uz = sinik * sinuk;
		double vx = xmx * cosuk - cosnok * sinuk;
		double vy = xmy * cosuk - sinnok * sinuk;
		double vz = sinik * cosuk;

		// Position
		double x = rk * ux;
		double y = rk * uy;
		double z = rk * uz;

		Vector vecPos(x, y, z);

		// Validate on altitude
		double altKm = (vecPos.Magnitude() * (XKMPER_WGS72 / AE));

		if (altKm < XKMPER_WGS72)
		{
			Julian decayTime = m_Orbit.Epoch();

			decayTime.AddMin(tsince);
			//throw cDecayException(decayTime, m_Orbit.SatName(true));
		}

		// Velocity
		double xdot = rdotk * ux + rfdotk * vx;
		double ydot = rdotk * uy + rfdotk * vy;
		double zdot = rdotk * uz + rfdotk * vz;

		Vector vecVel(xdot, ydot, zdot);

		Julian gmt = m_Orbit.Epoch();
		gmt.AddMin(tsince);

		m_Position = vecPos;
		m_Velocity = vecVel;
		m_Date = gmt;
	}

	// Orbital parameter variables which need only be calculated one
	// time for a given orbit (ECI position time-independent).
	double m_cosio;   double m_sinio;
	double m_betao2;  double m_betao;   double m_s4;
	double m_qoms24;  double m_tsi;     double m_eta;
	double m_eeta;    double m_coef;    double m_coef1;
	double m_c1;      double m_c3;      double m_c4;
	double m_a3ovk2;  double m_xmdot;   double m_omgdot;
	double m_xnodot;  double m_xnodcf;  double m_t2cof;

public:
	void ScalePosVector(double factor) { m_Position.Scale(factor); }
	void ScaleVelVector(double factor) { m_Velocity.Scale(factor); }

	const Vector& Position() const { return m_Position; }
	const Vector& Velocity() const { return m_Velocity; }
	Julian Date() const { return m_Date; }

public:
	const TLE& m_Orbit;

	Vector m_Position;
	Vector m_Velocity;
	Julian m_Date;
};
