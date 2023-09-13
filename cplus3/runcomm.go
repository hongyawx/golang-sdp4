package cplus3

/*
   #cgo CXXFLAGS: -std=c++11
   #include "comm.h"
*/
import "C"

func calculateAngle(line1, line2, line3 string, year, month, day, hour, min, sec int, lat, lon, altKm float64) AngleResult {
	cLine1 := C.CString(line1)
	cLine2 := C.CString(line2)
	cLine3 := C.CString(line3)

	/*defer func() {
		C.free(unsafe.Pointer(cLine1))
		C.free(unsafe.Pointer(cLine2))
		C.free(unsafe.Pointer(cLine3))
	}()*/

	result := C.calculateAngle(
		cLine1,
		cLine2,
		cLine3,
		C.int(year),
		C.int(month),
		C.int(day),
		C.int(hour),
		C.int(min),
		C.int(sec),
		C.double(lat),
		C.double(lon),
		C.double(altKm),
	)

	return AngleResult{
		Pitch: float32(result.pitch),
		Yaw:   float32(result.yaw),
	}
}

type AngleResult struct {
	Pitch float32
	Yaw   float32
}
