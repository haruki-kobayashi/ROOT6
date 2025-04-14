#ifndef VXXREADER_H
#define VXXREADER_H

#include <string>
#include <vector>
#include <array>
#include "netscan_data_types_ui.h"
#include "KeywordArgs.h"

#if defined VXX_DLL_BUILD
#define VXX_DLL_EXPORT __declspec(dllexport)
#else
#define VXX_DLL_EXPORT __declspec(dllimport)
#endif

namespace vxx
{

struct CutArea
{
	CutArea(double xmin, double xmax, double ymin, double ymax)
		: xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax),
		axmin(-DBL_MAX), axmax(DBL_MAX), aymin(-DBL_MAX), aymax(DBL_MAX)
	{};
	CutArea(double xmin, double xmax, double ymin, double ymax,
			double axmin, double axmax, double aymin, double aymax)
		: xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax),
		axmin(axmin), axmax(axmax), aymin(aymin), aymax(aymax)
	{};

	double	xmin;
	double	xmax;
	double	ymin;
	double	ymax;
	double	axmin;
	double	axmax;
	double	aymin;
	double	aymax;
};

struct HashEntry
{
	int32_t	index1;
	int32_t	index2;
	int	xmin;
	int	xmax;
	int	ymin;
	int	ymax;
	short	col_min;
	short	col_max;
	short	row_min;
	short	row_max;
};
struct ViewHeaderEntry
{
	ViewHeaderEntry() : RawId(-1) {};

	short	col;
	short	row;
	int	num;
	float	cx;
	float	cy;
	float	z1;
	float	z2;
	int32_t RawId;
};
struct GeometryInfo
{
	GeometryInfo()
		: pos(0), nseg(0), xmin(0), xmax(0), ymin(0), ymax(0),
		z(0), col_min(0), col_max(0), row_min(0), row_max(0),
		n(0), nPosInModule(0), nPlateInModule(0)
	{};

	int	pos;
	int	nseg;
	double	xmin, xmax;
	double	ymin, ymax;
	double	z;
	short	col_min, col_max;
	short	row_min, row_max;
	int	n;
	int	nPosInModule;
	int	nPlateInModule;
};

namespace opt
{

VXX_DEFINE_KEYWORD_OPTION_WITH_VALUE(index, VXX_TIE_ARGS(const std::array<int, 2>&))
VXX_DEFINE_KEYWORD_OPTION_WITH_VALUE(a, VXX_TIE_ARGS(const std::vector<CutArea>&))
VXX_DEFINE_KEYWORD_OPTION_WITH_VALUE(c, const std::string&)
VXX_DEFINE_KEYWORD_OPTION_WITH_VALUE(ph, int)

}

struct BvxxPolicy;
struct FvxxPolicy;

template <class TrackPolicy>
class VXX_DLL_EXPORT VxxReader
{
public:

	struct Parameters
	{
		Parameters(const std::string& filename, int pl, int zone)
			: mFilename(filename), mPlate(pl), mZone(zone), mPHMin(0), mIndex{ 0, -1 } {}
		std::string mFilename;
		int mPlate;
		int mZone;
		std::vector<CutArea> mAreaList;
		std::string mCorrmap;
		int mPHMin;
		std::array<int, 2> mIndex;
	};

	VxxReader()
		: mImpl(nullptr)
	{}
	~VxxReader();

	VxxReader(const VxxReader&) = delete;
	VxxReader(VxxReader&& v)
		: mImpl(v.mImpl)
	{
		v.mImpl = nullptr;
	}
	VxxReader& operator=(const VxxReader&) = delete;
	VxxReader& operator=(VxxReader&& v)
	{
		delete mImpl;
		mImpl = v.mImpl;
		v.mImpl = nullptr;
		return *this;
	}

	template <class ...Options, std::enable_if_t<AreAllKeywords<Options...>::value, std::nullptr_t> = nullptr>
	bool Begin(const std::string& filename, int pl, int zone, Options ...ops)
	{
		Parameters p(filename, pl, zone);
		SetOptions(p, ops...);
		return Begin(p);
	}
	bool Begin(const Parameters& p);
	bool NextHashEntry(HashEntry& h);
	template <bool B = std::is_same<TrackPolicy, BvxxPolicy>::value, std::enable_if_t<B, std::nullptr_t> = nullptr>
	bool NextBaseTrack(base_track_t& b);
	template <bool B = std::is_same<TrackPolicy, FvxxPolicy>::value, std::enable_if_t<B, std::nullptr_t> = nullptr>
	bool NextMicroTrack(micro_track_t& b);
	void End();

	template <class ...Options, std::enable_if_t<AreAllKeywords<Options...>::value, std::nullptr_t> = nullptr>
	decltype(auto) ReadAll(const std::string& filename, int pl, int zone, Options ...ops)
	{
		Parameters p(filename, pl, zone);
		SetOptions(p, ops...);
		return ReadAll(p);
	}
	template <bool B = std::is_same<TrackPolicy, BvxxPolicy>::value, std::enable_if_t<B, std::nullptr_t> = nullptr>
	std::vector<base_track_t> ReadAll(const Parameters& p);
	template <bool B = std::is_same<TrackPolicy, FvxxPolicy>::value, std::enable_if_t<B, std::nullptr_t> = nullptr, class Dummy = void>
	std::vector<micro_track_t> ReadAll(const Parameters& p);

	std::vector<ViewHeaderEntry> GetViewHeaderEntry(const std::string& filename, int pl, int zone);
	//Beginでファイルを開いた状態にある場合に限り、
	//こちらの関数を使ってもよい。
	std::vector<ViewHeaderEntry> GetViewHeaderEntry();

	GeometryInfo GetGeometryInfo(const std::string& filename, int pl, int zone);
	//Beginでファイルを開いた状態にある場合に限り、
	//こちらの関数を使ってもよい。
	GeometryInfo GetGeometryInfo();

private:

	class Impl;

	template <class ...Options, std::enable_if_t<AreAllKeywords<Options...>::value, std::nullptr_t> = nullptr>
	void SetOptions(Parameters& p, Options ...ops)
	{
		if constexpr (KeywordExists(opt::index, ops...)) p.mIndex = GetKeywordArg(opt::index, ops...);
		if constexpr (KeywordExists(opt::a, ops...)) p.mAreaList = GetKeywordArg(opt::a, ops...);
		if constexpr (KeywordExists(opt::c, ops...)) p.mCorrmap = GetKeywordArg(opt::c, ops...);
		if constexpr (KeywordExists(opt::ph, ops...)) p.mPHMin = GetKeywordArg(opt::ph, ops...);
	}
	template <class Type, bool B>
	struct NextTrack_impl;
	template <class Type, bool B>
	struct ReadAll_impl;

	Impl* mImpl;
};

namespace opt
{

VXX_DEFINE_KEYWORD_OPTION(append);
VXX_DEFINE_KEYWORD_OPTION_WITH_VALUE(views, const std::vector<ViewHeaderEntry>&);

}

template <class TrackPolicy>
class VXX_DLL_EXPORT VxxWriter
{
public:
	struct Parameters
	{
		Parameters(const std::string& filename, int pl, int zone)
			: mFilename(filename), mPlate(pl), mZone(zone), mAppend(false) {}
		std::string mFilename;
		int mPlate;
		int mZone;
		bool mAppend;
		std::vector<ViewHeaderEntry> mViews;
	};
	class Impl;

	VxxWriter()
		: mImpl(nullptr)
	{}
	~VxxWriter();

	VxxWriter(const VxxWriter&) = delete;
	VxxWriter(VxxWriter&& v)
		: mImpl(v.mImpl)
	{
		v.mImpl = nullptr;
	}
	VxxWriter& operator=(const VxxWriter&) = delete;
	VxxWriter& operator=(VxxWriter&& v)
	{
		delete mImpl;
		mImpl = v.mImpl;
		v.mImpl = nullptr;
		return *this;
	}

	template <class ...Options, std::enable_if_t<AreAllKeywords<Options...>::value, std::nullptr_t> = nullptr>
	bool Begin(const std::string& filename, int pl, int zone, Options ...ops)
	{
		Parameters p(filename, pl, zone);
		if constexpr (KeywordExists(opt::append, ops...)) p.mAppend = true;
		return Begin(p);
	}
	bool Begin(const Parameters& p);
	template <bool B = std::is_same<TrackPolicy, BvxxPolicy>::value, std::enable_if_t<B, std::nullptr_t> = nullptr>
	void Write(const base_track_t& b);
	template <bool B = std::is_same<TrackPolicy, FvxxPolicy>::value, std::enable_if_t<B, std::nullptr_t> = nullptr>
	void Write(const micro_track_t& b);
	void End();

	template <class Track, class ...Options, std::enable_if_t<AreAllKeywords<Options...>::value, std::nullptr_t> = nullptr>
	bool Write(const std::string& filename, int pl, int zone, const std::vector<Track>& t, Options ...ops)
	{
		Parameters p(filename, pl, zone);
		p.mAppend = GetKeywordArg(opt::append, ops..., false);
		if constexpr (std::is_same<Track, micro_track_t>::value)
			if constexpr (KeywordExists(opt::views, ops...))
				p.mViews = GetKeywordArg(opt::views, ops...);
			else
				p.mViews = CreateViews(t);
		else if constexpr (KeywordExists(opt::views, ops...))
			throw std::exception();//bvxxはViewHeaderEntryを必要としない。
		return Write(t, p);
	}
	template <bool B = std::is_same<TrackPolicy, BvxxPolicy>::value, std::enable_if_t<B, std::nullptr_t> = nullptr>
	bool Write(const std::vector<base_track_t>& t, const Parameters& p);
	template <bool B = std::is_same<TrackPolicy, FvxxPolicy>::value, std::enable_if_t<B, std::nullptr_t> = nullptr>
	bool Write(const std::vector<micro_track_t>& t, const Parameters& p);

	template <bool B = std::is_same<TrackPolicy, FvxxPolicy>::value, std::enable_if_t<B, std::nullptr_t> = nullptr>
	std::vector<ViewHeaderEntry> CreateViews(const std::vector<micro_track_t>& m) const;

private:

	template <class Type, bool B>
	struct Write_impl;

	Impl* mImpl;
};

using BvxxReader = VxxReader<BvxxPolicy>;
using FvxxReader = VxxReader<FvxxPolicy>;
using BvxxWriter = VxxWriter<BvxxPolicy>;
using FvxxWriter = VxxWriter<FvxxPolicy>;

}

#endif