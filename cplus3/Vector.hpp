//
// Vector.h
//
// Copyright 2003-2020 (c) Michael F. Henry
// Version 08/2020
//
#pragma once

#include "math.h"

class Vector
{
public:
	Vector(double x = 0.0, double y = 0.0, double z = 0.0, double w = 0.0)
		: m_x(x), m_y(y), m_z(z), m_w(w)
	{
	}

	virtual ~Vector() {};

	static double min(const double _Left, const double _Right)
	{
		return (_Right < _Left ? _Right : _Left);
	}

	static double max(const double& _Left, const double& _Right)
	{
		return (_Left < _Right ? _Right : _Left);
	}

	void Sub(const Vector& vec)
	{
		m_x -= vec.m_x;
		m_y -= vec.m_y;
		m_z -= vec.m_z;
		m_w -= vec.m_w;
	}

	void Scale(double factor)
	{
		m_x *= factor;
		m_y *= factor;
		m_z *= factor;
		m_w *= fabs(factor);
	}

	// angle between two vectors
	double Angle(const Vector& vec) const
	{
		double ratio = Dot(vec) / (Magnitude() * vec.Magnitude());

		// Avoid rounding errors
		if (ratio > 0.0) { ratio = min(ratio, 1.0f); }
		if (ratio < 0.0) { ratio = max(ratio, -1.0f); }

		return acos(ratio);
	}

	// vector magnitude
	double Magnitude() const
	{
		return sqrt((m_x * m_x) +
			(m_y * m_y) +
			(m_z * m_z));
	}

	// dot product
	double Dot(const Vector& vec) const
	{
		return (m_x * vec.m_x) +
			(m_y * vec.m_y) +
			(m_z * vec.m_z);
	}

	double m_x;
	double m_y;
	double m_z;
	double m_w;
};