#pragma once

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

		// Calculate the differences between the two basetracks
		// dax, day, dx, dy , dar, dal, dr, dl
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

	vxx::micro_track_subset_t ConvertNinjaToNetscan(const output_format_micro& micro);
	vxx::base_track_t ConvertNinjaToNetscan(const output_format_base& base);
	vxx::linklet_t ConvertNinjaToNetscan(const output_format_link& link);
	output_format_micro ConvertNetscanToNinja(const vxx::micro_track_subset_t& micro);
	output_format_base ConvertNetscanToNinja(const vxx::base_track_t& base);
	output_format_link ConvertNetscanToNinja(const vxx::linklet_t& link);
	bool PlausibleNetscan(const vxx::linklet_t& l);
	bool PlausibleNinja(const output_format_link& l);
}