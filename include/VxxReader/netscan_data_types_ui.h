#ifndef VXX_NETSCAN_DATA_TYPES_UI_H
#define VXX_NETSCAN_DATA_TYPES_UI_H

//
// This header-file can be used as an input to dump_linklet and l2c to select data members in binary-output. 
// To disable output, one can comment-out a member to be disabled. 
// Be careful that constructors will have to be corrected according to above changes or to be commented-out to supress compile-errors. 
//

#ifdef __CINT__
typedef Long64_t int64_t;
typedef Int_t int32_t;
typedef Short_t int16_t;
typedef UInt_t uint32_t;
typedef UShort_t uint16_t;
#else
#include <stdint.h>
#endif

namespace vxx {

//
// functions for alpha format
//
inline uint32_t hts_view_id( int col, int row ) { return ((uint32_t)(uint16_t)row&0xff00<<8)|((uint32_t)(uint16_t)col); };
inline uint32_t hts_imager_id( int /*col*/, int row ) { return ((uint32_t)(uint16_t)row&0x00ff); };

//
// functions for beta format
//
inline uint32_t hts_shot_id( int col, int row ) { return ((uint32_t)(uint16_t)row<<16)|((uint32_t)(uint16_t)col); };

class micro_track_t {
public:
	micro_track_t() : pos(0),col(0),row(0),zone(0),isg(0),rawid(0),ph(0),ax(0),ay(0),x(0),y(0),z(0),z1(0),z2(0),px(0),py(0) {};
	micro_track_t( int pos, int col, int row, int zone, int isg , int64_t rawid, int ph, double ax, double ay, double x, double y, double z, double z1, double z2, float px, float py ) 
		: pos(pos),col(col),row(row),zone(zone),isg(isg),rawid(rawid),ph(ph),ax(ax),ay(ay),x(x),y(y),z(z),z1(z1),z2(z2),px(px),py(py) {};

	double ax,ay;
	double x,y,z;
	double z1,z2;
	float px,py;
	int ph,pos,col,row,zone,isg;
	int64_t rawid;
	uint64_t gid()const {
		return ((uint64_t)(uint16_t)row) << 48 | ((uint64_t)(uint16_t)col) << 32 | (uint64_t)isg;
	}
};

class micro_track_subset_t {
public:
	micro_track_subset_t() : ph(0),ax(0),ay(0),z(0),pos(0),col(0),row(0),zone(0),isg(0),rawid(0) {};
	micro_track_subset_t( int ph, double ax, double ay, double z, int pos, int col, int row, int zone, int isg, int64_t rawid ) 
		: ph(ph),ax(ax),ay(ay),z(z),pos(pos),col(col),row(row),zone(zone),isg(isg),rawid(rawid) {};

	double ax,ay;
	double z;
	int ph;
	int pos,col,row,zone,isg;
	int64_t rawid;
	uint64_t gid()const {
		return ((uint64_t)(uint16_t)row) << 48 | ((uint64_t)(uint16_t)col) << 32 | (uint64_t)isg;
	}
};

class base_track_t {
public:
	base_track_t() : pl(0),rawid(0),isg(0),zone(0),ax(0),ay(0),x(0),y(0),z(0) {};
	base_track_t ( int pl, int64_t rawid, int isg, int zone, double ax, double ay, double x, double y, double z, const micro_track_subset_t& m1, const micro_track_subset_t& m2 ) 
		: pl(pl),rawid(rawid),isg(isg),zone(zone),ax(ax),ay(ay),x(x),y(y),z(z) { m[0]=m1; m[1]=m2; };
    base_track_t(double x, double y, double ax, double ay, double z_begin, double z_end, int phvol) {
        this->ax = ax;
        this->ay = ay;
        m[0].ax = ax;
        m[0].ay = ay;
        m[1].ax = ax;
        m[1].ay = ay;
        this->x = x;
        this->y = y;
        m[0].ph = phvol;
        m[1].ph = phvol;
        z = (z_begin + z_end) * 0.5;
        m[0].z = z_begin;
        m[1].z = z_end;
    }
	double ax,ay;
	double x,y,z;
	int pl;
	int isg,zone;
	int dmy;		// dummy member to fix bug in root interpreter mode ( 4byte align for member access, while sizeof(base_track_t) returns 176 i.e. 8byte align )
	int64_t rawid;
	std::array<micro_track_subset_t, 2> m;
	double dangle()const {
		return std::sqrt(std::pow((m[0].ax - ax), 2) + std::pow((m[0].ay - ay), 2) + std::pow((m[1].ax - ax), 2) + std::pow((m[1].ay - ay), 2));
	}
	double cx()const {
		return  x + ax * (m[1].z - m[0].z) * 0.5;
	}
	double cy()const {
		return  y + ay * (m[1].z - m[0].z) * 0.5;
	}
	double x1()const {
		return  x + ax * (m[1].z - m[0].z) * 1.0;
	}
	double y1()const {
		return  y + ay * (m[1].z - m[0].z) * 1.0;
	}
	double x0()const {
		return  x;
	}
	double y0()const {
		return  y;
	}
	void trans(double aff_coef[6]) {
		double bx = aff_coef[0] * x + aff_coef[1] * y + aff_coef[4];
		double by = aff_coef[2] * x + aff_coef[3] * y + aff_coef[5];
		double bax = aff_coef[0] * ax + aff_coef[1] * ay;
		double bay = aff_coef[2] * ax + aff_coef[3] * ay;
		x = bx;
		y = by;
		ax = bax;
		ay = bay;
		for (int i = 0; i < 2; i++) {
			bax = aff_coef[0] * m[i].ax + aff_coef[1] * m[i].ay;
			bay = aff_coef[2] * m[i].ax + aff_coef[3] * m[i].ay;
			m[i].ax = bax;
			m[i].ay = bay;
		}
	}
};

class linklet_t {
public:
	linklet_t() : zproj(0),dx(0),dy(0),xc(0),yc(0) { pos[0] = pos[1] = 0; };

	int pos[2];
	double zproj;
	double dx,dy;
	double xc,yc;
	std::array<base_track_t, 2> b;
};

class sub_base_track {
public:
    int _cx, _cy; //0.01um
    int _ax, _ay; //0.01mrad
    int _dangle; //0.01mrad
    int phvol;
    uint64_t gid;
    void set_cxcy(double cx, double cy) {
        _cx = (int)(cx * 100);
        _cy = (int)(cy * 100);
    }
    double cx() const { return _cx * 0.01; }
    double cy() const { return _cy * 0.01; }
    double ax() const { return _ax * 0.00001; }
    double ay() const { return _ay * 0.00001; }
    double dangle() const { return _dangle * 0.00001; }
    sub_base_track() {}
    sub_base_track(int cx, int cy, int ax, int ay, int dangle, int phvol, uint64_t gid) :
        _cx(cx), _cy(cy), _ax(ax), _ay(ay), _dangle(dangle), phvol(phvol), gid(gid) {};
    void trans(double aff_coef[6]) {
        double bx = aff_coef[0] * cx() + aff_coef[1] * cy() + aff_coef[4];
        double by = aff_coef[2] * cx() + aff_coef[3] * cy() + aff_coef[5];
        double bax = aff_coef[0] * ax() + aff_coef[1] * ay();
        double bay = aff_coef[2] * ax() + aff_coef[3] * ay();
        _cx = (int)(bx * 100);
        _cy = (int)(by * 100);
        _ax = (int)(bax * 100000);
        _ay = (int)(bay * 100 * 1000);
    }
};


};

#endif
