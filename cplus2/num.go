package cplus2

// #cgo LDFLAGS: -L. -lstdc++
// #cgo CXXFLAGS: -std=c++14 -I.
// #include "num.h"
import "C"
import "unsafe"

/**
LDFLAGS: -L. -lstdc++ 表示 在当前的目录中搜索目录查找库文件 链接器在链接库时使用 libstdc++库
#cgo CXXFLAGS: -std=c++14 -I.  表示在当前目录查找头文件
*/

type GoNum struct {
	num C.Num
}

func New() GoNum {
	var ret GoNum
	ret.num = C.NumInit()
	return ret
}
func (n GoNum) Free() {
	C.NumFree((C.Num)(unsafe.Pointer(n.num)))
}
func (n GoNum) Inc() {
	C.NumIncrement((C.Num)(unsafe.Pointer(n.num)))
}
func (n GoNum) GetValue() int {
	return int(C.NumGetValue((C.Num)(unsafe.Pointer(n.num))))
}
