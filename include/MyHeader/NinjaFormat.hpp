#pragma once

#include <cmath>

#include <VxxReader/netscan_data_types_ui.h>
#include <MyHeader/CalcTrackPair.hpp>

namespace NinjaFormat
{
	class output_format_micro {
	public:
		int pos, view, imager, zone, isg, ph, vph, px;
	};

	class output_format_base {
	public:
		int pl, rawid;
		double ax, ay, x, y, z;
		output_format_micro m[2];
	};

	class output_format_link {
	public:
		output_format_base b[2];
		double dax, day, dx, dy, dar, dal, dr, dl;

		void Calc_difference() {
			dax = b[0].ax - b[1].ax;
			day = b[0].ay - b[1].ay;
			dx = b[0].x - b[1].x;
			dy = b[0].y - b[1].y;
			double tan_theta = std::sqrt(b[0].ax*b[0].ax + b[0].ay*b[0].ay);
			dar = (dax*b[0].ax + day*b[0].ay) / tan_theta;
			dal = (dax*b[0].ay - day*b[0].ax) / tan_theta;
			dr = CalcTrackPair::RadialDistance(b[0], b[1]);
			dl = CalcTrackPair::LateralDistance(b[0], b[1]);
		};
	};

	vxx::micro_track_subset_t ConvertNinjaToNetscan(const output_format_micro& micro) {
		vxx::micro_track_subset_t m;

		constexpr uint8_t NumberOfImager = 72;
		uint32_t ShotID = micro.view * NumberOfImager + micro.imager;

		m.pos = micro.pos;
		m.col = (int16_t)(ShotID&0x0000ffff);
		m.row = (int16_t)((ShotID&0xffff0000)>>16);
		m.zone = micro.zone;
		m.isg = micro.isg;
		m.rawid = 0; // not used
		m.ph = micro.ph * 10000 + micro.vph;
		m.ax = 0.0; // not used
		m.ay = 0.0; // not used
		m.z = 0.0; // not used

		return m;
	}

	vxx::base_track_t ConvertNinjaToNetscan(const output_format_base& base) {
		vxx::base_track_t b;

		b.pl = base.pl;
		b.rawid = base.rawid;
		b.zone = 0; // not used
		b.isg = 0; // not used
		b.ax = base.ax;
		b.ay = base.ay;
		b.x = base.x;
		b.y = base.y;
		b.z = base.z;

		b.m[0] = ConvertNinjaToNetscan(base.m[0]);
		b.m[1] = ConvertNinjaToNetscan(base.m[1]);

		return b;
	}

	vxx::linklet_t ConvertNinjaToNetscan(const output_format_link& link) {
		vxx::linklet_t l;

		l.pos[0] = link.b[0].pl * 10;
		l.pos[1] = link.b[1].pl * 10;
		l.zproj = l.xc = l.yc = 0.;
		l.dx = link.dx;
		l.dy = link.dy;

		l.b[0] = ConvertNinjaToNetscan(link.b[0]);
		l.b[1] = ConvertNinjaToNetscan(link.b[1]);

		return l;
	}

	output_format_micro ConvertNetscanToNinja(const vxx::micro_track_subset_t& micro) {
		output_format_micro m;

		constexpr uint8_t NumberOfImager = 72;
		uint32_t ShotID = ((uint32_t)(uint16_t)micro.row << 16) | ((uint32_t)(uint16_t)micro.col);

		m.view = ShotID / NumberOfImager;
		m.imager = ShotID % NumberOfImager;

		m.pos = micro.pos;
		m.zone = micro.zone;
		m.isg = micro.isg;
		m.ph = micro.ph / 10000;
		m.vph = micro.ph % 10000;
		m.px = 0; // not used

		return m;
	}

	output_format_base ConvertNetscanToNinja(const vxx::base_track_t& base) {
		output_format_base b;

		b.pl = base.pl;
		b.rawid = base.rawid;
		b.ax = base.ax;
		b.ay = base.ay;
		b.x = base.x;
		b.y = base.y;
		b.z = base.z;

		b.m[0] = ConvertNetscanToNinja(base.m[0]);
		b.m[1] = ConvertNetscanToNinja(base.m[1]);

		return b;
	}

	output_format_link ConvertNetscanToNinja(const vxx::linklet_t& link) {
		output_format_link l;

		l.b[0] = ConvertNetscanToNinja(link.b[0]);
		l.b[1] = ConvertNetscanToNinja(link.b[1]);

		l.Calc_difference();

		return l;
	}

	bool PlausibleNetscan(const vxx::linklet_t& l) {
		if (l.b[0].m[0].ph > 16 || l.b[0].m[1].ph > 16) return false;
		if (l.b[1].m[0].ph > 16 || l.b[1].m[1].ph > 16) return false;
		if (l.b[0].m[0].ph < 1 || l.b[0].m[1].ph < 1) return false;
		if (l.b[1].m[0].ph < 1 || l.b[1].m[1].ph < 1) return false;
		if (l.pos[0] < 0 || l.pos[0] > 10000) return false;
		if (l.pos[1] < 0 || l.pos[1] > 10000) return false;
		if (std::abs(l.b[0].ax) > 100. || std::abs(l.b[0].ay) > 100.) return false;
		if (std::abs(l.b[1].ax) > 100. || std::abs(l.b[1].ay) > 100.) return false;
		if (std::abs(l.b[0].z) > 1e7 || std::abs(l.b[1].z) > 1e7) return false;
		if (std::abs(l.dx) > 1e7 || std::abs(l.dy) > 1e7) return false;
		return true;
	}

	bool PlausibleNinja(const output_format_link& l) {
		if (l.b[0].m[0].ph > 16 || l.b[0].m[1].ph > 16) return false;
		if (l.b[1].m[0].ph > 16 || l.b[1].m[1].ph > 16) return false;
		if (l.b[0].m[0].ph < 1 || l.b[0].m[1].ph < 1) return false;
		if (l.b[1].m[0].ph < 1 || l.b[1].m[1].ph < 1) return false;
		if (l.b[0].pl < 0 || l.b[0].pl > 1000) return false;
		if (l.b[1].pl < 0 || l.b[1].pl > 1000) return false;
		if (std::abs(l.b[0].ax) > 100. || std::abs(l.b[0].ay) > 100.) return false;
		if (std::abs(l.b[1].ax) > 100. || std::abs(l.b[1].ay) > 100.) return false;
		if (std::abs(l.b[0].z) > 1e7 || std::abs(l.b[1].z) > 1e7) return false;
		if (std::abs(l.dx) > 1e7 || std::abs(l.dy) > 1e7) return false;
		return true;
	}
}