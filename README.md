# Golang calls c++ to calculate elevation and azimuth using the sdp4 algorithm

## Description

## Installation



## Example

```go

func TestSGp4(t *testing.T) {
// 输入 tle 3行格式数据 + 时间 ，计算一个 仰角 和方位
result := calculateAngle("STARLINK-1011", "1 44717U 19074E   23251.77703217  .00022635  00000+0  15319-2 0  9999", "2 44717  53.0543 344.8465 0001241  87.7302 272.3829 15.06440356211049", 2023, 9, 13, 21, 27, 23, 30.6, 104.03, 0.5)
fmt.Println(result)
//输出结果为 {-21.096296 294.50558}
}

```