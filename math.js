/* 
* 常用空间分析函数
* 包含 点、线、面 之间的相互关系
* 本文主要记录一些自己开发中常用的一些分析函数，比较正宗和全面的可以看一些空间分析库，比如前端的 Turf 和 JSTS
*/

// 点 -------------------------------------------------------------------------

// 点到点的距离
function dist2d(coord1, coord2) {
  let dx = coord1[0] - coord2[0];
  let dy = coord1[1] - coord2[1];
  return Math.sqrt(dx * dx + dy * dy)
}
dist2d([1, 1], [4, 5]) // 5


//判断两个点集是否相等
function equals(coord1, coord2) {
  let equals = true;
  for (let i = coord1.length - 1; i >= 0; --i){
    if (coord1[i] != coord2[i]) {
      equals = false;
      break
    }
  }
  return equals
}
equals([1, 3], [1, 3]) // true
equals([1, 3, 1], [1, 3, 2]) // false


// 线 -------------------------------------------------------------------------

// 线段的长度
function formatLength(coords) {
  coords = coords || []
  let length = 0
  //通过遍历坐标计算两点之前距离，进而得到整条线的长度
  for (let i = 0, leng = coords.length - 1; i < leng; i++) {
    length += dist2d(coords[i], coords[i + 1]);
  }
  return length
}
formatLength([[1, 0], [2, 1], [3, 0], [4, 2]]) // 5.06449510224598


// 根据偏移量偏移一条线段
function lineOffset(coords, deltaX, deltaY){
  deltaX = deltaX || 0
  deltaX = isNaN(deltaX) ? 0 : deltaX
  deltaY = deltaY || 0
  deltaY = isNaN(deltaY) ? 0 : deltaY

  if(deltaX == 0 && deltaY == 0)return coords

  coords.forEach(coord => {
    coord[0] += deltaX;
    coord[1] += deltaY;
  })
  return coords
}
lineOffset([[1, 1], [2, 3], [3, 3], [4, 2], [2, 0]], 2, 3)
//         [[3, 4], [4, 6], [5, 6], [6, 5], [4, 3]]


// 线段上距离点P最近的一个点
function getNearestCoord(point, lines) {
  let d, res = { dist: Infinity }
	if(!(Array.isArray(lines[0]) && Array.isArray(lines[0][0]))){
		lines = [lines]
	}
	lines.forEach(function(coords){
		for (let i = 0; i < coords.length; i++) {
			d = dist2d(point, coords[i]);
			if (d < res.dist) {
        res.dist = d;
        res.index = i
        res.point = coords[i]
			}
		}
	})
	
	return res.point ? res : null;
}
getNearestCoord([2.2, 3.1], [[1, 1], [2, 3], [3, 3], [4, 2], [2, 0]])
// {dist: 0.5099019513592785, index: 1, point: [2, 3]}

/* 
* 点P到线段AB的最短距离
* 使用矢量算法，计算线AP在线段AB方向上的投影，当需要计算的数据量很大时，这种方式优势明显
* 特殊情况如点在线段上、点在端点、点在线段延长线上等等的情况全部适用于此公式，只是作为特殊情况出现，无需另作讨论。这也是矢量算法思想的优势所在。
* 函数返回值：point 投影坐标  dist 点P到投影距离  type 垂足位置，不为0表示垂足在线段外
*/
function pointToSegmentDist(point, point1, point2){
  let x = point[0], x1 = point1[0], x2 = point2[0]
  let y = point[1], y1 = point1[1], y2 = point2[1]

  //线段AB 为一个点
  if(x1 == x2 && y1 == y2) return {
    type: 0,
    point: point1,
    dist: 0
  }

  let cross = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1);
  //let r = cross / d2
  //r < 0 点P的垂足在线段AB外，且点P距离线段AB最近的点为A
  //r = 0 点P的垂足和点P距离线段AB最近的点为A
  if (cross <= 0) return {
    type: 1,
    point: point1,
    dist: Math.sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1))
  };
    
  let d2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
  //r > 1 点P的垂足在线段AB外，且点P距离线段AB最近的点为B
  //r = 1 点P的垂足和点P距离线段AB最近的点为B
  if (cross >= d2) return {
    type: 2,
    point: point2,
    dist: Math.sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2))
  };
    
  let r = cross / d2;
  let px = x1 + (x2 - x1) * r;
  let py = y1 + (y2 - y1) * r;
  return {
    type: 0,
    point: [px, py],
    dist: Math.sqrt((x - px) * (x - px) + (py - y) * (py - y))
  };
}
pointToSegmentDist([1, 3], [0, 1], [3, 1]); 
//{ dist: 2, point: [1, 1], type: 0 }

/* 
* pointToSegmentDist变种：点P到直线AB的垂足及最短距离
*/
function pointToFootDist(point, point1, point2){
  let x = point[0], x1 = point1[0], x2 = point2[0]
  let y = point[1], y1 = point1[1], y2 = point2[1]

  //线段AB 为一个点
  if(x1 == x2 && y1 == y2) return {
    type: true,
    point: point1,
    dist: 0
  }

  let dx = x2 - x1;
  let dy = y2 - y1;
  //r < 0 点P的垂足在线段AB外，且点P距离线段AB最近的点为A
  //r > 1 点P的垂足在线段AB外，且点P距离线段AB最近的点为B
  //r = 0 点P的垂足和点P距离线段AB最近的点为A
  //r = 1 点P的垂足和点P距离线段AB最近的点为B
  let r = (dx * (x0 - x1) + dy * (y0 - y1)) / (dx * dx + dy * dy || 0);

  let px = x1 + dx * r;
  let py = y1 + dy * r;
  return {
    type: r >= 0 && r <= 1, //true  垂足在线段内   false 垂足在线段外
    point: [px, py], //垂足
    dist: Math.sqrt((x - px) * (x - px) + (py - y) * (py - y)) //点P到垂足距离
  };
}
pointToFootDist([3, 1], [1, 0], [3, 0])
// {dist: 1, point: [3, 0], type: true}


/*
* 计算点P到直线AB的距离
* 使用三角形面积，计算线AP在直线AB方向上的投影
* Area   = |(1/2)(x1 * y2 + x2 * y3 + x3 * y1 - x2 * y1 - x3 * y2 - x1 * y3)|
* Bottom = Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))
* Area = .5 * Bottom * H
* Height = Area / .5 / Bottom
*/
function verticalDistance(point, point1, point2){
  if(equals(point1, point2)){
    return dist2d(point, point1)
  }
  let area = Math.abs(0.5 * (point1[0] * point2[1] + point2[0] * point[1] + point[0] * point1[1] - point2[0] * point1[1] - point[0] * point2[1] - point1[0] * point[1]));
  let bottom = dist2d(point1, point2)
  let height = area / 0.5 / bottom;

  return height;
}
verticalDistance([4, 2], [1, 0], [3, 0]) // 2

/* 
* 计算点P到直线AB的距离
* 根据公式，计算线AP在直线AB方向上的投影
* d = Math.abs((A * x + B * y + C) / Math.sqrt( A * A + B * B))
* 公式介绍：https://zhuanlan.zhihu.com/p/26307123
* 直线公式：y = kx + b
* 直线公式：kx - y + b = 0
* 假设 B = -1，则 A = k，C = b   
* 直线斜率：k = (y1 - y2) / (x1 - x2)
* 常数计算：b = y - kx
*/
function verticalDistance2(point, point1, point2){
	//如果point1[0] == point2[0] 说明是条竖着的线
	if(point1[0] == point2[0]){
		return Math.abs(point[0] - point1[0])
	}else{
    let k = (point1[1] - point2[1]) / (point1[0] - point2[0])
    let b = point1[1] - k * point1[0]
		return Math.abs((k * point[0] - point[1] + b) / Math.sqrt( k * k + 1))
	}
}
verticalDistance2([4, 2], [1, 0], [3, 0]) // 2


/* 
* 计算纬度和经度指定的两点之间的距离(米)
* 在计算涉及到经纬度的计算时(长度、面积等)，如果只是粗略的计算，那么我们使用坐标点计算，将计算结果按照每经纬度111km转换即可。这个结果越往两极失真越大
* 如果想要精细点，就需要用到地球长半轴参数去计算
* WGS84 坐标系下的长半轴参数 radius = 6378137 m
* WGS84(World Geodetic System 1984) 坐标系是为GPS全球定位系统使用而建立的坐标系统。是国际公认的坐标系
* 北京54坐标系下的长半轴参数 radius = 6378245 m
* 1954年北京坐标系可以认为是前苏联1942年坐标系的延伸。它的原点不在北京而是在前苏联的普尔科沃。所以误差较大，缺点较多。
* 西安80坐标系下的长半轴参数 radius = 6378140 ± 5 m
* 1980西安坐标系所采用的IAG1975椭球，其长半轴要比WGS84椭球长半轴的值大3米左右，而这可能引起地表长度误差达10倍左右
* CGCS2000(2000 国家大地坐标系) 长半轴参数 radius = 6378137 m
* 2000国家大地坐标系，是我国当前最新的国家大地坐标系，是全球地心坐标系在我国的具体体现，其原点为包括海洋和大气的整个地球的质量中心。
* 正弦曲线投影下的长半轴参数 radius = 6370997 m
* 正弦曲线投影是一种等面积的伪圆柱投影。规定纬线投影为平行直线，经线投影为对称于中央经线的正弦曲线，同一纬线上经距相等，纬距向两极缩小。主要用于小比例尺世界地图
*/
function haversineDistance(c1, c2) {
  let radius = 6378137
  let lat1 = toRadians(c1[1]);
  let lat2 = toRadians(c2[1]);
  let deltaLatBy2 = (lat2 - lat1) / 2;
  let deltaLonBy2 = toRadians(c2[0] - c1[0]) / 2;
  let a = Math.sin(deltaLatBy2) * Math.sin(deltaLatBy2) + Math.sin(deltaLonBy2) * Math.sin(deltaLonBy2) * Math.cos(lat1) * Math.cos(lat2);
  return 2 * radius * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a))
}
function toRadians(_){
  return _ * Math.PI / 180
}
haversineDistance([117, 34], [118, 39])  //563732.8911197125米
haversineDistance([118, 39], [117, 34])  //563732.8911197125米


/* 
* 根据一个点切割线段
* 点不必落在线段上
* 常用于计算鼠标位置距离线段起点或终点的距离
* 简单版的是使用 getNearestCoord，而不是 pointToSegmentDist
* getNearestCoord 寻找的是线上真实存在的一点
* pointToSegmentDist 寻找的点位于线上，但不一定真实存在与线上
*/
function getClosestPoint(coords, point){
  if(!coords || coords.length < 2) return [[], []]

  let squaredDistance = {dist: Infinity}, index
  for(let i = 1; i < coords.length; i++){
    let d = pointToSegmentDist(point, coords[i - 1], coords[i])
    if(d.dist < squaredDistance.dist){
      squaredDistance = d
      index = i
    }
  }

  if(index === undefined){
    return [coords, []]
  }

  let prearr = coords.slice(0, index)
  if(prearr.length && !equals(squaredDistance.point, prearr[prearr.length - 1])){
    prearr.push(squaredDistance.point)
  }

  let nextarr = coords.slice(index)
  if(nextarr.length && !equals(squaredDistance.point, nextarr[0])){
    nextarr.unshift(squaredDistance.point)
  }

  return [prearr, nextarr]
}
getClosestPoint([[0, 1], [2, 3], [5, 6], [3, 4], [4, 4]], [3, 5])
/*
[
	[[0,1],[2,3],[3.5,4.5]],
	[[3.5,4.5],[5,6],[3,4],[4,4]]
]
*/


/* 
* 线段与线段的交点
* 解线性方程组, 求线段AB与线段CD的交点坐标，如果没有交点，返回null
*/
function intersects(coords1, coords2) {
  let x1 = coords1[0][0];
  let y1 = coords1[0][1];
  let x2 = coords1[1][0];
  let y2 = coords1[1][1];
  let x3 = coords2[0][0];
  let y3 = coords2[0][1];
  let x4 = coords2[1][0];
  let y4 = coords2[1][1];
  //斜率交叉相乘 k1 = (y4 - y3) / (x4 - x3)    k2 = (y2 - y1) / (x2 - x1)
  //k1 k2 同乘 (x4 - x3) * (x2 - x1) 并相减得到denom
  let denom = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1)); 
  // 如果分母为0 则平行或共线, 不相交 
  if (denom === 0) {
    return null;
  }

  let numeA = ((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3));
  let numeB = ((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3));
  let uA = numeA / denom;
  let uB = numeB / denom;

  // 交点在线段1上，且交点也在线段2上  
  if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
    let x = x1 + (uA * (x2 - x1));
    let y = y1 + (uA * (y2 - y1));
    return [x, y];
  }
  return null;
}
intersects([[0,0],[1,1]], [[3,0],[2,1]])  //null
intersects([[0,0],[1,1]], [[3,0],[0,1]])  //[0.75, 0.75]


/* 
* 线与线的所有交点
* 接收两条线的点集合，求取line2 和 line1的交点
* 如果规定count参数，代表返回前多少个交点，否则全部返回
* 返回交点集合，集合中的元素属性包括 index1：交点在line1的位置，index2：交点在line2的位置，coords：交点坐标
* 代码参考truf.lineIntersect：http://turfjs.org/docs/#lineIntersect
* 如果交点为线的点，则会重复返回，包括返回的数据结构，都是为下面切割面函数splitPolygon服务的
*/
function lineIntersect(line1, line2, count){
  let result = []
	for(let i = 1; i < line1.length; i++){
    let coords1 = [line1[i-1], line1[i]]
    //求取数据的边界范围，函数放在了下面的多边形中
    let bbox1 = bbox(coords1)
		for(let j = 1; j < line2.length; j++){
      let coords2 = [line2[j-1], line2[j]]
      let bbox2 = bbox(coords2)
      //判断两个边界范围的关系： bbox1 是否包含 bbox2，函数放在了下面的多边形中
      if(isIntersects(bbox1, bbox2)){
        let p = intersects(coords1, coords2)
        if(p){
          result.push({index1: i, index2: j, coords: p})
          if(count && (result.length >= count)){
            return result
          }
        }
      }
		}
  }
	return result
}
lineIntersect([[0, 0], [0, 1], [1, 3], [3, 2], [5, 0]], [[-1, 0], [4, 3]])
//[{"index1":1,"index2":1,"coords":[0,0.6]},{"index1":3,"index2":1,"coords":[2.6363636363636367,2.1818181818181817]}]


/* 
* 判断线段路径是顺时针还是逆时针
* 使用格林公式计算
* 返回值 d < 0  '顺时针' : d > 0  '逆时针'
* truf对应实现truf.booleanClockwise://turfjs.org/docs/#booleanClockwise
*/
function isClockwise(line){
  let length = line.length
  let d = 0

  if(length.length < 3)return d;

  //如果不是闭合线，则改为闭合线
  if(!equals(line[0], line[length - 1])){
    length = (line = line.slice()).push(line[0]);
  }

  //沿着多边形的边求曲线积分,若积分为正,则是沿着边界曲线正方向(逆时针),反之为顺时针
  //最后一个点是开始点，只参与末尾点(倒数第二个点)的计算，不单独计算
  for(let i = 0; i < length - 1; i++){
    d += -0.5 * (line[i + 1][0] - line[i][0]) * (line[i + 1][1] + line[i][1])
  }

  return d
}
isClockwise([[0,0],[1,1],[1,0],[0,0]]) // -0.5 顺时针
isClockwise([[0,0],[1,0],[1,1],[0,0]]) //  0.5 逆时针

/* 
* 判断线段路径是顺时针还是逆时针
* 使用端点判断
* 为正时，p1-p2-p3   路径的走向为逆时针，  
* 为负时，p1-p2-p3   走向为顺时针，  
* 为零时，p1-p2-p3   所走的方向不变，亦即三点在一直线上
* 返回值 d < 0  '顺时针' : d > 0  '逆时针'
*/
function isClockwise2(line){
  let length = line.length
  let d = 0

  if(length.length < 3)return d;

  //如果不是闭合线，则改为闭合线
  if(!equals(line[0], line[length - 1])){
    length = (line = line.slice()).push(line[0]);
  }

  //循环遍历多边形的坐标选取X或者Y值中最大或者最小的点，这个点必然是凸点
  //然后取该点前后各一个点Pm-1、Pm+1，组成向量(Pm-1,Pm)、(Pm,Pm+1)。然后进行向量叉乘即可判断出顺时针或逆时针
  let max = line[0][0], maxIndex = 0;
  for(let i = 1; i < length - 1; i++){
    if(line[i][0] > max){
      maxIndex = i
      max = line[i][0]
    }
  }

  //找到端点 p2 ，然后取改点前后各一点 p1 p3
  let p1 = line[maxIndex ? (maxIndex - 1) : (length - 2)]
  let p2 = line[maxIndex]
  let p3 = line[maxIndex + 1]
  //然后根据三个点组成的两个向量乘积判断，进行向量叉乘即可判断出顺时针或逆时针
  //d = (x2 - x1) * (y3 - y2) - (x3 - x2) * (y2 - y1)
  d = (p2[0] - p1[0]) * (p3[1] - p2[1]) - (p3[0] - p2[0]) * (p2[1] - p1[1])

  return d
}
isClockwise2([[0,0],[1,1],[1,0],[0,0]]) // -1 顺时针
isClockwise2([[0,0],[1,0],[1,1],[0,0]]) //  1 逆时针


/* 
* 简化线 与 平滑线
* 简化线 与 平滑线 比较复杂，这里就不展示代码了
* 简化线 又可以叫做线的抽稀，一般使用的方法为 道格拉斯-普克，这也是我使用的方法
* 道格拉斯-普克算法(Douglas–Peucker algorithm，亦称为拉默-道格拉斯-普克算法、迭代适应点算法、分裂与合并算法)是将曲线近似表示为一系列点，并减少点的数量的一种算法。其思想主要是保留关键点。

* Douglas-Peucker算法描述：
* 1、在线首尾两点A，B之间连接一条直线AB，该直线为曲线的弦；
* 2、得到曲线上离该直线段距离最大的点C，计算其与AB的距离d；
* 3、比较该距离与预先给定的阈值tolerance的大小，如果小于tolerance，则该直线段作为曲线的近似，该段曲线处理完毕。
* 4、如果距离大于阈值tolerance，则用C将曲线分为两段AC和BC，并分别对两段曲线进行1~3的处理。
* 5、当所有曲线都处理完毕时，依次连接各个分割点形成的折线，即可以作为曲线的近似。

* 除了可以使用道格拉斯-普克算法对线进行抽稀外，还可以使用Wang-Müller算法(保留关键折弯)、Zhou-Jones算法(保留加权有效面积)、Visvalingam-Whyatt算法(保留有效面积)等。
* 在实际应用中发现，现实是残酷的，如果阈值过小，简化后的控制点会很多，用户很难进行操作。而阈值变大，控制点变少了，但是经过平滑之后，线面变形严重，比较失真。真实和美观之间很难平衡，难以兼得。


* 平滑线 一般使用的方法为 贝赛尔曲线插值和B样条曲线插值
* 在应用中发现，两者在某些时候都会有些问题，所以自己对贝赛尔曲线算法做了一点改变，具体的可以看https://segmentfault.com/a/1190000031626358
*/



// 多边形 -------------------------------------------------------------------------

// 获取多边形边界范围
function bbox(coords) {
  // x/经度最小值 y/纬度最小值 x/经度最大值 y/纬度最大值
  let res = [Infinity, Infinity, -Infinity, -Infinity];
  coords.forEach(coord => {
    if (res[0] > coord[0]) res[0] = coord[0];
    if (res[2] < coord[0]) res[2] = coord[0];
    if (res[1] > coord[1]) res[1] = coord[1];
    if (res[3] < coord[1]) res[3] = coord[1];
  })
  return res;
}
bbox([[1, 1], [2, 3], [3, 3], [4, 2], [2, 0]])
// [1, 0, 4, 3]


//判断两个边界范围的关系： a 是否包含 b
function isContains(a, b) {
  return a[0] <= b[0] &&
    a[1] <= b[1] &&
    b[2] <= a[2] &&
    b[3] <= a[3];
}
isContains([1, 0, 4, 3], [2, 2, 5, 5])  //false
isContains([1, 0, 4, 3], [2, 1, 3, 2])  //true


//判断两个边界范围的关系： a 与 b 是否有交集
function isIntersects(a, b) {
  return b[0] <= a[2] &&
    b[1] <= a[3] &&
    b[2] >= a[0] &&
    b[3] >= a[1];
}
isIntersects([1, 0, 4, 3], [2, 2, 5, 5])  //true
isIntersects([1, 0, 4, 3], [5, 2, 5, 5])  //false


// 获取多边形中心
function center(coords) {
  let ext = bbox(coords);
  let x = (ext[0] + ext[2]) / 2;
  let y = (ext[1] + ext[3]) / 2;
  return [x, y];
}
center([[1, 1], [2, 3], [3, 3], [4, 2], [2, 0]])
// [2.5, 1.5]


// 获取多边形重心
function centroid(coords) {
  let xSum = 0, ySum = 0, len = 0;
  coords.forEach(coord => {
    xSum += coord[0];
    ySum += coord[1];
    len++;
  })

  return [xSum / len, ySum / len];
}
centroid([[1, 1], [2, 3], [3, 3], [4, 2], [2, 0]])
// [2.4, 1.8]


// 获取多边形质心
function centerOfMass(coords) {
  let centr = centroid(coords);
  let neutralizedPoints = coords.map(function (point$$1) {
    return [
      point$$1[0] - centr[0],
      point$$1[1] - centr[1]
    ];
  });
  
  let sx = 0, sy = 0, sArea = 0;
  let x1, x2, y1, y2, a;
  for (let i = 0, len = coords.length - 1; i < len; i++) {
    x1 = neutralizedPoints[i][0]; 
    y1 = neutralizedPoints[i][1];
    x2 = neutralizedPoints[i + 1][0]; 
    y2 = neutralizedPoints[i + 1][1];

    // a 是计算有符号面积和最终坐标的公因子
    a = x1 * y2 - x2 * y1;
    // sArea 用来计算有符号面积的和
    sArea += a;
    // sx 和 sy 是用于计算最终坐标的总和
    sx += (x1 + x2) * a;
    sy += (y1 + y2) * a;
  }

  if (sArea === 0) {
    return centre;
  } else {
    // 计算有符号面积，并因式分解: x = 1 / 6A
    let area = sArea * 0.5;
    let areaFactor = 1 / (6 * area);

    return [
      centr[0] + areaFactor * sx,
      centr[1] + areaFactor * sy
    ]
  }
}
centerOfMass([[1, 1], [2, 3], [3, 3], [4, 2], [2, 0]])
// [2.569230769230769, 1.8735042735042735]


/* 
* 获取多边形面积
* 原理：△ABC的面积就是'向量AB'和'向量AC'两个向量叉积的绝对值的一半。其正负表示三角形顶点是在右手系还是左手系。
* 所以可以把多边形以多边形内一点为顶点，拆分出一个个三角形，再求取面积
* 设 n 边形的顶点依次是 (x1，y1)(x2,y2)......(xn,yn)
* 那么：s = (x1y2-x2y1)/2 + (x2y3-x3y2)/2 +......+ (xny1-x1yn)/2
* 获取的 s 是有向面积，需要取绝对值
*/
function getArea(polygon){
  let area = 0;
  let len = polygon.length
  let x1 = polygon[len - 1][0]
  let y1 = polygon[len - 1][1]
  for(let i = 0; i < len; i++){
    let x2 = polygon[i][0];
    let y2 = polygon[i][1];
    area += y1 * x2 - x1 * y2;
    x1 = x2;
    y1 = y2
  }
  return Math.abs(area / 2)
}
getArea([[0,0],[0,10],[10,10],[10,0]])  //100
getArea([[0,0],[0,3],[4,0]])  //6


/* 
* 计算多边形投影到地球上时的近似面积（平方米）
* 注意，如果环是顺时针方向的，这个区域将是正的，否则它将是负的
* 注意，投影及其长半轴参数radius的介绍这里不做叙述，可以看上面haversineDistance的相关介绍
*/
function geodesicArea(coords) {
  let radius = 6378137
  let area = 0, len = coords.length;
  if(len <= 2){ return area }

  let x1 = coords[len - 1][0];
  let y1 = coords[len - 1][1];
  for (let i = 0; i < len; i++) {
    let x2 = coords[i][0], y2 = coords[i][1];
    area += toRadians(x2 - x1) * (2 + Math.sin(toRadians(y1)) + Math.sin(toRadians(y2)));
    x1 = x2;
    y1 = y2
  }
  area = area * radius * radius / 2
  return Math.abs(area)
}
geodesicArea([[118, 39], [117, 34], [117, 33], [116, 36], [117, 40]])  //64519860945.307144平方米


/* 
* 判断点是否在多边形内
* 射线法（ray casting）或者奇偶规则法（even odd rule）
* 从这个点做一条射线，计算它跟多边形边界的交点个数，如果交点个数为奇数，那么点在多边形内部，否则点在多边形外。
* 对于射线法，需要排除以下几种情况：
* 1、点在多边形的边上
* 2、点和多边形的顶点重合
* 3、射线经过多边形顶点
* 4、射线刚好经过多边形的一条边
* truf对应实现truf.booleanPointInPolygon://turfjs.org/docs/#booleanPointInPolygon
*/
function pointInPolygon(point, polygon) {
  let px = point[0], py = point[1], flag = false

  let sx, sy, tx, ty
  for(let i = 0, l = polygon.length, j = l - 1; i < l; j = i, i++) {
    sx = polygon[i][0]
    sy = polygon[i][1]
    tx = polygon[j][0]
    ty = polygon[j][1]

    // 点与多边形顶点重合       
    if((sx === px && sy === py) || (tx === px && ty === py)) {
      return true
    }

    // 判断线段两端点是否在射线两侧       
    if((sy < py && ty >= py) || (sy >= py && ty < py)) {
      // 线段上与射线 Y 坐标相同的点的 X 坐标         
      let x = sx + (py - sy) * (tx - sx) / (ty - sy)

      // 点在多边形的边上         
      if(x === px) {
        return true
      }

      // 射线穿过多边形的边界         
      if(x > px) {
        flag = !flag
      }
    }
  }

  // 射线穿过多边形边界的次数为奇数时点在多边形内     
  return flag
}

/* 
* 判断点是否在多边形内
* 回转数法
* 用线段分别连接点和多边形的全部顶点，计算所有点与相邻顶点连线的夹角，计算所有夹角和，最后根据角度累加值计算回转数
* 注意每个夹角都是有方向的，所以有可能是负值。360°（2π）相当于一次回转。
* 当回转数为 0 时，点在闭合曲线外部。
*/
function pointInPolygon2(point, polygon) {
  let px = point[0], py = point[1], sum = 0

  let sx, sy, tx, ty
  for(let i = 0, l = polygon.length, j = l - 1; i < l; j = i, i++) {
    sx = polygon[i][0]
    sy = polygon[i][1]
    tx = polygon[j][0]
    ty = polygon[j][1]

    // 点与多边形顶点重合       
    if((sx === px && sy === py) || (tx === px && ty === py)) {
      return true
    }

    // 点与相邻顶点连线的夹角       
    let angle = Math.atan2(sy - py, sx - px) - Math.atan2(ty - py, tx - px)

    // 确保夹角不超出取值范围（-π 到 π）
    if(angle >= Math.PI) {
      angle = angle - Math.PI * 2
    } else if(angle <= -Math.PI) {
      angle = angle + Math.PI * 2
    }

    sum += angle
  }

  // 计算回转数并判断点和多边形的几何关系     
  return Math.round(sum / Math.PI) !== 0
}

/* 
* 判断点是否在多边形内
* openlayers中的方法
*/
function pointInPolygon3(point, polygon) {
  let x = point[0], y = point[1]
  let wn = 0;
  let x1 = polygon[0][0];
  let y1 = polygon[0][1];
  for (let i = 0; i < polygon.length; i++) {
    let x2 = polygon[i][0];
    let y2 = polygon[i][1];
    if (y1 <= y) {
      if (y2 > y && (x2 - x1) * (y - y1) - (x - x1) * (y2 - y1) > 0){
        wn++;
      }
    } else if (y2 <= y && (x2 - x1) * (y - y1) - (x - x1) * (y2 - y1) < 0){
      wn--;
    }
    x1 = x2;
    y1 = y2
  }
  return wn !== 0
}


/* 
* 两个多边形之间的关系：相交、包含、相离
* https://segmentfault.com/a/1190000020916225
* 1、遍历A多边形上的点，判断是否有坐标点在B多边形内 --- 返回结果 a
* 2、遍历B多边形上的点，判断是否有坐标点在A多边形内 --- 返回结果 b
* 如果a、b都为true，则两个多边形相交
* 如果a为true，b为false，则多边形B包含多边形A
* 如果a为false，b为true，则多边形A包含多边形B
* 如果a、b都为false，则两个多边形远离
* truf对应实现 包含 truf.booleanContains://turfjs.org/docs/#booleanContains
* truf对应实现 交叉 truf.booleanCrosses://turfjs.org/docs/#booleanCrosses
* truf对应实现 相离 truf.booleanDisjoint://turfjs.org/docs/#booleanDisjoint
*/
function judge(coordsA, coordsB){
  let boola = coordsA.some(item => {
    return pointInPolygon3(item, coordsB)
  }) ? 1 : 0;
  let boolb = coordsB.some(item => {
    return pointInPolygon3(item, coordsA)
  }) ? 1 : 0;

  return ['相离', 'A包含B', 'B包含A', '相交'][boola * 2 + boolb]
}


/* 
* 线切割多边形(闭合线)
* 使用线将多边形切割为一个个小多边形
* truf中类似的有：
* 差异 通过从第一个多边形中剪切第二个多边形来找到两个多边形之间的差异。
* truf.difference://turfjs.org/docs/#difference
* 相交 取两个多边形，并找到它们的交点。如果它们相交，则返回相交边界。如果它们不相交，则返回undefined。
* truf.intersect://turfjs.org/docs/#intersect
*/
function splitPolygonByLine(polygon, line){
  let polygonItem = { coords: polygon }
  _splitPolygonByLine(polygonItem, line)
  return getChildrenPolygons([polygonItem])
}
function _splitPolygonByLine(polygonItem, lineCoords){
  if(polygonItem.children && polygonItem.children.length){
    return polygonItem.children.forEach(polygon => {
      _splitPolygonByLine(polygon, lineCoords)
    });
  }

  lineCoords = lineCoords.slice();
  let points = lineIntersect(lineCoords, polygonItem.coords, 2)
  //console.log('points', points)
  if(points.length >= 2){
    //判断线的开始点在面内还是面外，以判断哪两个点练成的线在面外，舍弃在外面的
    let startOut = pointInPolygon3(lineCoords[0], polygonItem.coords)
    let point1Out = pointInPolygon3(lineCoords[points[0].index1], polygonItem.coords)
    if(startOut && !point1Out){
      lineCoords = lineCoords.slice(points[0].index1);
      return _splitPolygonByLine(polygonItem, lineCoords)
    }

    let newPolygons = splitPolygon(polygonItem, lineCoords, points);
    //console.log('newPolygons', newPolygons)
    if(newPolygons.length >= 2){
      polygonItem.children = newPolygons

      lineCoords = lineCoords.slice(points[1].index1)
      newPolygons.forEach(polygon => {
        _splitPolygonByLine(polygon, lineCoords)
      })
    }
  }
}

function splitPolygon(polygonItem, lineCoords, points){
  let result = [];

  //lineCoords = lineCoords.slice()
  let polygonCoords = polygonItem.coords
  
  let startIndex = 0, endIndex = 0, lineIndex = 0

  points.sort((a,b)=>a.index1 - b.index1)
  
  points.reduce((prePoint, point, index) => {
    startIndex = prePoint.index1
    endIndex = point.index1

    let line;
    if(endIndex == startIndex){
      //return point
      line = [prePoint.coords, point.coords]
    }else{
      line = lineCoords.slice(startIndex, endIndex)
      line.unshift(prePoint.coords)
      line.push(point.coords)
    };

    lineIndex = endIndex

    let p1, p2 = polygonCoords.slice();
    if(prePoint.index2 > point.index2){
      startIndex = point.index2
      endIndex = prePoint.index2

      p1 = polygonCoords.slice(startIndex, endIndex)
      p1 = p1.concat(line)
      p1.push(p1[0])

      line.reverse()
      p2 = p2.slice(0, startIndex).concat(line).concat( p2.slice(endIndex) )
    }else if(prePoint.index2 < point.index2){
      startIndex = prePoint.index2
      endIndex = point.index2

      p2 = p2.slice(0, startIndex).concat(line).concat( p2.slice(endIndex) )

      line.reverse()
      p1 = polygonCoords.slice(startIndex, endIndex)
      p1 = p1.concat(line)
      p1.push(p1[0])
    }else{
      //startIndex = endIndex = prePoint.index1
      p1 = line.slice()

      let pc = polygonItem.coords[point.index2 - 1]
      if(dist2d(prePoint.coords, pc) > dist2d(point.coords, pc)){
        line.reverse()
      }

      p2.splice(point.index2, 0, ...line)
    }

    result.push({ coords: p1 }, { coords: p2 })

    return point
  })

  return result
}
function getChildrenPolygons(polygons){
  let result = []
  polygons.forEach(item => {
    if(item.children){
      result = result.concat( getChildrenPolygons(item.children) )
    }else if(getArea(item.coords) > 0){
      result.push(item.coords)
    }
  });
  
  return result
}
splitPolygonByLine([[0, 0], [0, 10], [10, 10], [10, 0], [0, 0]], [[-1, 0], [6, 12], [12, 2], [-1, 9]])
/* 
[
  [[0,1.7142857142857142],[0,8.461538461538462],[2.995121951219513,6.848780487804878],[0,1.7142857142857142]],
  [[0,10],[4.833333333333334,10],[2.995121951219513,6.848780487804878],[0,8.461538461538462],[0,10]],
  [[10,10],[10,5.333333333333334],[7.2,10],[10,10]],
  [[4.833333333333334,10],[7.2,10],[10,5.333333333333334],[10,3.076923076923077],[2.9951219512195113,6.848780487804879],[4.833333333333334,10]],
  [[0,0],[0,1.7142857142857142],[2.9951219512195113,6.848780487804879],[10,3.076923076923077],[10,0],[0,0]]
]
*/
