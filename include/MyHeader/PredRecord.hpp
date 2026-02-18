#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <string>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <map>
#include <unordered_map>
#include <optional>
#include <cmath>

namespace pred {

/* ===============================
   data structures
   =============================== */

struct Pred {
    int64_t chain_id;
    int pl;
    double ax, ay;
    double x, y, z;
    int ph;
    double mean_ax, mean_ay;
    double dal_sigma;
    int pl0, pl1;
    int nseg;
    /*
        ax, ay, x, y    : 隣接プレート（2 or 1枚）のbasetrackから計算される角度・位置
        mean_ax, mean_ay: chain内basetrackの平均角度
        dal_sigma       : chain内basetrackの平均角度に対する各basetrackのlateral角度ずれの標準偏差
                          小さいほど直線的なchainであることを示す
        pl0, pl1        : chain内の最初と最後のプレート番号
        nseg            : chain内のセグメント数
    */
};

struct Reco {
    int64_t rawid; // basetrack rawid
    int pl;
    double ax, ay;
    double x, y, z;
    int ph;
};

struct Record {
    Pred pred;
    std::vector<Reco> recos; // 複数 reco を許容
};

// 一行分のデータ
struct RowData {
    Pred pred;
    std::optional<Reco> reco;
};

/* ===============================
   derived quantities
   =============================== */

inline bool found(const Record& r) {
    return !r.recos.empty();
};

inline int nreco(const Record& r) {
    return static_cast<int>(r.recos.size());
};

inline bool found(const RowData& r) {
    return r.reco.has_value();
};

/* ===============================
   column schema
   =============================== */

enum class ColumnType { Int64, Int, Double, Bool };

using WriteFn = void (*)(std::ostream&, const Record&, const Reco*);
using ReadFn = void (*)(const std::string&, RowData&, Reco&, bool&);
using GetInt64Fn = int64_t (*)(const Record&, const Reco*);
using GetIntFn = int (*)(const Record&, const Reco*);
using GetDoubleFn = double (*)(const Record&, const Reco*);
using GetBoolFn = bool (*)(const Record&, const Reco*);
using SetInt64Fn = void (*)(int64_t, RowData&, Reco&, bool&);
using SetIntFn = void (*)(int, RowData&, Reco&, bool&);
using SetDoubleFn = void (*)(double, RowData&, Reco&, bool&);
using SetBoolFn = void (*)(bool, RowData&, Reco&);

struct ColumnSchema {
    const char* name;
    ColumnType  type;
    WriteFn write;
    ReadFn  read;
    GetInt64Fn get_int64;
    GetIntFn get_int;
    GetDoubleFn get_double;
    GetBoolFn get_bool;
    SetInt64Fn set_int64;
    SetIntFn set_int;
    SetDoubleFn set_double;
    SetBoolFn set_bool;
};

inline const char* to_string(ColumnType t) {
    switch (t) {
        case ColumnType::Int64: return "int64_t";
        case ColumnType::Int:    return "int";
        case ColumnType::Double: return "double";
        case ColumnType::Bool:   return "bool";
    }
    return "unknown";
}

constexpr int64_t MISSING_INT64  = -INT64_MAX;
constexpr int      MISSING_INT  = -INT_MAX;
inline double missing_double() {
    return std::numeric_limits<double>::quiet_NaN();
}

inline void write_missing(std::ostream& o, ColumnType t) {
    switch (t) {
        case ColumnType::Int64: o << MISSING_INT64;      break;
        case ColumnType::Int:    o << MISSING_INT;      break;
        case ColumnType::Double: o << missing_double(); break;
        case ColumnType::Bool:   o << false;            break;
    }
}

struct StreamFormatGuard {
    std::ostream& os;
    std::ios::fmtflags flags;
    std::streamsize precision;

    StreamFormatGuard(std::ostream& o)
        : os(o), flags(o.flags()), precision(o.precision()) {}

    ~StreamFormatGuard() {
        os.flags(flags);
        os.precision(precision);
    }

    StreamFormatGuard(const StreamFormatGuard&) = delete;
    StreamFormatGuard& operator=(const StreamFormatGuard&) = delete;
};

inline void write_fixed_precision(std::ostream& o, double v, int prec) {
    StreamFormatGuard guard(o);
    o << std::fixed << std::setprecision(prec) << v;
}

// Binary format constants
constexpr uint32_t BINARY_MAGIC = 0x50524543;  // "PREC"
constexpr uint16_t BINARY_VERSION = 1;

/* ===============================
   schema table (single source of truth)
   =============================== */

inline constexpr std::array<ColumnSchema, 24> SCHEMA = {{
    {"pred.chain_id", ColumnType::Int64,
     [](auto& o, auto& r, auto*) { o << r.pred.chain_id; },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.chain_id = std::stoull(token); },
     [](const Record& r, const Reco*) { return r.pred.chain_id; },
     nullptr, nullptr, nullptr,
     [](int64_t v, RowData& row, Reco&, bool&) { row.pred.chain_id = v; },
     nullptr, nullptr, nullptr},

    {"pred.pl", ColumnType::Int,
     [](auto& o, auto& r, auto*) { o << r.pred.pl; },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.pl = std::stoi(token); },
     nullptr,
     [](const Record& r, const Reco*) { return r.pred.pl; },
     nullptr, nullptr, nullptr,
     [](int v, RowData& row, Reco&, bool&) { row.pred.pl = v; },
     nullptr, nullptr},

    {"pred.ax", ColumnType::Double,
     [](auto& o, auto& r, auto*) { write_fixed_precision(o, r.pred.ax, 4); },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.ax = std::stod(token); },
     nullptr, nullptr,
     [](const Record& r, const Reco*) { return r.pred.ax; },
     nullptr, nullptr, nullptr,
     [](double v, RowData& row, Reco&, bool&) { row.pred.ax = v; },
     nullptr},

    {"pred.ay", ColumnType::Double,
     [](auto& o, auto& r, auto*) { write_fixed_precision(o, r.pred.ay, 4); },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.ay = std::stod(token); },
     nullptr, nullptr,
     [](const Record& r, const Reco*) { return r.pred.ay; },
     nullptr, nullptr, nullptr,
     [](double v, RowData& row, Reco&, bool&) { row.pred.ay = v; },
     nullptr},

    {"pred.x",  ColumnType::Double,
     [](auto& o, auto& r, auto*) { write_fixed_precision(o, r.pred.x, 1); },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.x = std::stod(token); },
     nullptr, nullptr,
     [](const Record& r, const Reco*) { return r.pred.x; },
     nullptr, nullptr, nullptr,
     [](double v, RowData& row, Reco&, bool&) { row.pred.x = v; },
     nullptr},

    {"pred.y",  ColumnType::Double,
     [](auto& o, auto& r, auto*) { write_fixed_precision(o, r.pred.y, 1); },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.y = std::stod(token); },
     nullptr, nullptr,
     [](const Record& r, const Reco*) { return r.pred.y; },
     nullptr, nullptr, nullptr,
     [](double v, RowData& row, Reco&, bool&) { row.pred.y = v; },
     nullptr},

    {"pred.z",  ColumnType::Double,
     [](auto& o, auto& r, auto*) { write_fixed_precision(o, r.pred.z, 1); },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.z = std::stod(token); },
     nullptr, nullptr,
     [](const Record& r, const Reco*) { return r.pred.z; },
     nullptr, nullptr, nullptr,
     [](double v, RowData& row, Reco&, bool&) { row.pred.z = v; },
     nullptr},

    {"pred.ph", ColumnType::Int,
     [](auto& o, auto& r, auto*) { o << r.pred.ph; },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.ph = std::stoi(token); },
     nullptr,
     [](const Record& r, const Reco*) { return r.pred.ph; },
     nullptr, nullptr, nullptr,
     [](int v, RowData& row, Reco&, bool&) { row.pred.ph = v; },
     nullptr, nullptr},

    {"pred.mean_ax", ColumnType::Double,
     [](auto& o, auto& r, auto*) { write_fixed_precision(o, r.pred.mean_ax, 4); },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.mean_ax = std::stod(token); },
     nullptr, nullptr,
     [](const Record& r, const Reco*) { return r.pred.mean_ax; },
     nullptr, nullptr, nullptr,
     [](double v, RowData& row, Reco&, bool&) { row.pred.mean_ax = v; },
     nullptr},

    {"pred.mean_ay", ColumnType::Double,
     [](auto& o, auto& r, auto*) { write_fixed_precision(o, r.pred.mean_ay, 4); },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.mean_ay = std::stod(token); },
     nullptr, nullptr,
     [](const Record& r, const Reco*) { return r.pred.mean_ay; },
     nullptr, nullptr, nullptr,
     [](double v, RowData& row, Reco&, bool&) { row.pred.mean_ay = v; },
     nullptr},

    {"pred.dal_sigma", ColumnType::Double,
     [](auto& o, auto& r, auto*) { write_fixed_precision(o, r.pred.dal_sigma, 6); },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.dal_sigma = std::stod(token); },
     nullptr, nullptr,
     [](const Record& r, const Reco*) { return r.pred.dal_sigma; },
     nullptr, nullptr, nullptr,
     [](double v, RowData& row, Reco&, bool&) { row.pred.dal_sigma = v; },
     nullptr},

    {"pred.pl0", ColumnType::Int,
     [](auto& o, auto& r, auto*) { o << r.pred.pl0; },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.pl0 = std::stoi(token); },
     nullptr,
     [](const Record& r, const Reco*) { return r.pred.pl0; },
     nullptr, nullptr, nullptr,
     [](int v, RowData& row, Reco&, bool&) { row.pred.pl0 = v; },
     nullptr, nullptr},

    {"pred.pl1", ColumnType::Int,
     [](auto& o, auto& r, auto*) { o << r.pred.pl1; },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.pl1 = std::stoi(token); },
     nullptr,
     [](const Record& r, const Reco*) { return r.pred.pl1; },
     nullptr, nullptr, nullptr,
     [](int v, RowData& row, Reco&, bool&) { row.pred.pl1 = v; },
     nullptr, nullptr},

    {"pred.nseg", ColumnType::Int,
     [](auto& o, auto& r, auto*) { o << r.pred.nseg; },
     [](const auto& token, auto& row, auto&, auto&) { row.pred.nseg = std::stoi(token); },
     nullptr,
     [](const Record& r, const Reco*) { return r.pred.nseg; },
     nullptr, nullptr, nullptr,
     [](int v, RowData& row, Reco&, bool&) { row.pred.nseg = v; },
     nullptr, nullptr},

    // derived quantity (policy-generated)
    {"pred.found", ColumnType::Bool,
     [](auto& o, auto& r, auto*) { o << found(r); },
     [](const auto&, auto&, auto&, auto&) { /* derived, skip */ },
     nullptr, nullptr, nullptr,
     [](const Record& r, const Reco*) { return found(r); },
     nullptr, nullptr, nullptr, nullptr},
    // number of reco tracks associated to this pred
    {"pred.nreco", ColumnType::Int,
     [](auto& o, auto& r, auto*) { o << nreco(r); },
     [](const auto&, auto&, auto&, auto&) { /* derived, skip */ },
     nullptr,
     [](const Record& r, const Reco*) { return nreco(r); },
     nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},

    {"reco.rawid", ColumnType::Int64,
     [](auto& o, auto&, auto* rc) {
        if (rc) o << rc->rawid;
        else write_missing(o, ColumnType::Int64);
     },
     [](const auto& token, auto&, auto& reco, auto& has_reco) {
        reco.rawid = (token == std::string(std::to_string(MISSING_INT64)))
            ? MISSING_INT64 : std::stoll(token);
        if (reco.rawid != MISSING_INT64) has_reco = true;
     },
     [](const Record&, const Reco* rc) { return rc ? rc->rawid : MISSING_INT64; },
     nullptr, nullptr, nullptr,
     [](int64_t v, RowData&, Reco& reco, bool& has_reco) {
        reco.rawid = v;
        if (v != MISSING_INT64) has_reco = true;
     },
     nullptr, nullptr, nullptr},

    {"reco.pl", ColumnType::Int,
     [](auto& o, auto&, auto* rc) {
        if (rc) o << rc->pl;
        else write_missing(o, ColumnType::Int);
     },
     [](const auto& token, auto&, auto& reco, auto& has_reco) {
        reco.pl = (token == std::to_string(MISSING_INT))
            ? MISSING_INT : std::stoi(token);
        if (reco.pl != MISSING_INT) has_reco = true;
     },
     nullptr,
     [](const Record&, const Reco* rc) { return rc ? rc->pl : MISSING_INT; },
     nullptr, nullptr, nullptr,
     [](int v, RowData&, Reco& reco, bool& has_reco) {
        reco.pl = v;
        if (v != MISSING_INT) has_reco = true;
     },
     nullptr, nullptr},

    {"reco.ax", ColumnType::Double,
     [](auto& o, auto&, auto* rc) {
        if (rc) write_fixed_precision(o, rc->ax, 4);
        else write_missing(o, ColumnType::Double);
     },
     [](const auto& token, auto&, auto& reco, auto& has_reco) {
        reco.ax = (token == "nan") ? missing_double() : std::stod(token);
        if (!std::isnan(reco.ax)) has_reco = true;
     },
     nullptr, nullptr,
     [](const Record&, const Reco* rc) { return rc ? rc->ax : missing_double(); },
     nullptr, nullptr, nullptr,
     [](double v, RowData&, Reco& reco, bool& has_reco) {
        reco.ax = v;
        if (!std::isnan(v)) has_reco = true;
     },
     nullptr},

    {"reco.ay", ColumnType::Double,
     [](auto& o, auto&, auto* rc) {
        if (rc) write_fixed_precision(o, rc->ay, 4);
        else write_missing(o, ColumnType::Double);
     },
     [](const auto& token, auto&, auto& reco, auto& has_reco) {
        reco.ay = (token == "nan") ? missing_double() : std::stod(token);
        if (!std::isnan(reco.ay)) has_reco = true;
     },
     nullptr, nullptr,
     [](const Record&, const Reco* rc) { return rc ? rc->ay : missing_double(); },
     nullptr, nullptr, nullptr,
     [](double v, RowData&, Reco& reco, bool& has_reco) {
        reco.ay = v;
        if (!std::isnan(v)) has_reco = true;
     },
     nullptr},

    {"reco.x", ColumnType::Double,
     [](auto& o, auto&, auto* rc) {
        if (rc) write_fixed_precision(o, rc->x, 1);
        else write_missing(o, ColumnType::Double);
     },
     [](const auto& token, auto&, auto& reco, auto& has_reco) {
        reco.x = (token == "nan") ? missing_double() : std::stod(token);
        if (!std::isnan(reco.x)) has_reco = true;
     },
     nullptr, nullptr,
     [](const Record&, const Reco* rc) { return rc ? rc->x : missing_double(); },
     nullptr, nullptr, nullptr,
     [](double v, RowData&, Reco& reco, bool& has_reco) {
        reco.x = v;
        if (!std::isnan(v)) has_reco = true;
     },
     nullptr},

    {"reco.y", ColumnType::Double,
     [](auto& o, auto&, auto* rc) {
        if (rc) write_fixed_precision(o, rc->y, 1);
        else write_missing(o, ColumnType::Double);
     },
     [](const auto& token, auto&, auto& reco, auto& has_reco) {
        reco.y = (token == "nan") ? missing_double() : std::stod(token);
        if (!std::isnan(reco.y)) has_reco = true;
     },
     nullptr, nullptr,
     [](const Record&, const Reco* rc) { return rc ? rc->y : missing_double(); },
     nullptr, nullptr, nullptr,
     [](double v, RowData&, Reco& reco, bool& has_reco) {
        reco.y = v;
        if (!std::isnan(v)) has_reco = true;
     },
     nullptr},

    {"reco.z", ColumnType::Double,
     [](auto& o, auto&, auto* rc) {
        if (rc) write_fixed_precision(o, rc->z, 1);
        else write_missing(o, ColumnType::Double);
     },
     [](const auto& token, auto&, auto& reco, auto& has_reco) {
        reco.z = (token == "nan") ? missing_double() : std::stod(token);
        if (!std::isnan(reco.z)) has_reco = true;
     },
     nullptr, nullptr,
     [](const Record&, const Reco* rc) { return rc ? rc->z : missing_double(); },
     nullptr, nullptr, nullptr,
     [](double v, RowData&, Reco& reco, bool& has_reco) {
        reco.z = v;
        if (!std::isnan(v)) has_reco = true;
     },
     nullptr},

    {"reco.ph", ColumnType::Int,
     [](auto& o, auto&, auto* rc) {
        if (rc) o << rc->ph;
        else write_missing(o, ColumnType::Int);
     },
     [](const auto& token, auto&, auto& reco, auto& has_reco) {
        reco.ph = (token == std::to_string(MISSING_INT))
            ? MISSING_INT : std::stoi(token);
        if (reco.ph != MISSING_INT) has_reco = true;
     },
     nullptr,
     [](const Record&, const Reco* rc) { return rc ? rc->ph : MISSING_INT; },
     nullptr, nullptr, nullptr,
     [](int v, RowData&, Reco& reco, bool& has_reco) {
        reco.ph = v;
        if (v != MISSING_INT) has_reco = true;
     },
     nullptr, nullptr},
}};

/* ===============================
   table writer
   =============================== */

class FileWriter {
public:
    std::vector<const ColumnSchema*> columns;
    std::string delimiter;
    std::string format;  // "csv", "tsv", "binary"
    bool format_explicit = false;  // 明示的にフォーマットが指定されたか

    // デフォルトコンストラクタ: 全カラムを出力（CSVモード）
    FileWriter() : delimiter(","), format("csv"), format_explicit(false) {
        for (const auto& s : SCHEMA) {
            columns.push_back(&s);
        }
    }

    // コンストラクタ: フォーマットを指定（"csv", "tsv", または "binary"）
    FileWriter(const std::string& fmt)
        : delimiter(","), format("csv"), format_explicit(true) {
        if (fmt == "tsv" || fmt == "\\t") {
            delimiter = "\t";
            format = "tsv";
        } else if (fmt == "binary" || fmt == "bin") {
            format = "binary";
        } else if (fmt == "csv" || fmt == ",") {
            // デフォルトのまま
        } else {
            throw std::runtime_error(
                "Invalid format specified: " + fmt);
        }
        for (const auto& s : SCHEMA) {
            columns.push_back(&s);
        }
    }

    // コンストラクタ: カスタムカラムを指定して初期化（CSVモード）
    FileWriter(const std::vector<std::string>& column_names)
        : delimiter(","), format("csv"), format_explicit(false) {
        SetColumns(column_names);
    }

    // コンストラクタ: カスタムカラムとフォーマットを指定
    FileWriter(const std::vector<std::string>& column_names, const std::string& fmt)
        : delimiter(","), format("csv"), format_explicit(true) {
        if (fmt == "tsv" || fmt == "\\t") {
            delimiter = "\t";
            format = "tsv";
        } else if (fmt == "binary" || fmt == "bin") {
            format = "binary";
        } else if (fmt == "csv" || fmt == ",") {
            // デフォルトのまま
        } else {
            throw std::runtime_error(
                "Invalid format specified: " + fmt);
        }
        SetColumns(column_names);
    }

    // フォーマットを設定
    void SetFormat(const std::string& fmt) {
        if (fmt == "tsv" || fmt == "\\t") {
            delimiter = "\t";
            format = "tsv";
        } else if (fmt == "binary" || fmt == "bin") {
            format = "binary";
        } else if (fmt == "csv" || fmt == ",") {
            delimiter = ",";
            format = "csv";
        } else {
            throw std::runtime_error(
                "Invalid format specified: " + fmt);
        }
        format_explicit = true;  // 明示的な指定をマーク
    }

    // カスタムカラムを設定する
    void SetColumns(const std::vector<std::string>& column_names) {
        columns.clear();
        for (const auto& name : column_names) {
            const ColumnSchema* s = find_schema(name);
            if (!s) {
                throw std::runtime_error(
                    "Unknown column name: " + name);
            }
            columns.push_back(s);
        }
    }

    void WriteHeader(std::ostream& out) const {
        if (!out) {
            throw std::runtime_error("Output stream is not valid");
        }
        if (columns.empty()) {
            throw std::runtime_error("No columns defined for output");
        }

        if (format == "binary") {
            // バイナリヘッダー: マジックナンバーとバージョン
            uint32_t magic = BINARY_MAGIC;
            uint16_t version = BINARY_VERSION;
            out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
            out.write(reinterpret_cast<const char*>(&version), sizeof(version));

            // カラム情報
            uint32_t num_cols = static_cast<uint32_t>(columns.size());
            out.write(reinterpret_cast<const char*>(&num_cols), sizeof(num_cols));

            for (const auto& col : columns) {
                uint32_t name_len = static_cast<uint32_t>(std::strlen(col->name));
                out.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
                out.write(col->name, name_len);
                uint8_t type = static_cast<uint8_t>(col->type);
                out.write(reinterpret_cast<const char*>(&type), sizeof(type));
            }
        } else {
            // テキストヘッダー
            out << "# format: pred-reco-table v1\n";
            out << "# missing:\n";
            out << "#   int64_t = " << std::to_string(MISSING_INT64) << "\n";
            out << "#   int     = " << std::to_string(MISSING_INT) << "\n";
            out << "#   double  = nan\n";
            out << "# columns:\n# ";

            for (size_t i = 0; i < columns.size(); ++i) {
                out << columns[i]->name << "|" << to_string(columns[i]->type);
                if (i + 1 < columns.size()) out << delimiter;
            }
            out << "\n";
        }
    }

    // バイナリ形式のときはすべてのRecordをまとめて書き込む必要がある
    // テキスト形式のときは任意の単位でを書き込める
    void WriteRecords(std::ostream& out, const std::vector<Record>& records) const {
        if (!out) {
            throw std::runtime_error("Output stream is not valid");
        }
        if (columns.empty()) {
            throw std::runtime_error("No columns defined for output");
        }

        // 同じ (chain_id, pl) を持つ Record をまとめる
        std::map<std::pair<int64_t, int>, Record> map;
        for (const auto& r : records) {
            auto key = std::make_pair(r.pred.chain_id, r.pred.pl);
            auto& rec = map[key];
            if (rec.recos.empty()) {
                rec.pred = r.pred;  // 初回のみ pred を設定
            }
            // すべての reco を追加
            for (const auto& rc : r.recos) {
                rec.recos.push_back(rc);
            }
        }

        if (format == "binary") {
            // バイナリ形式
            // pred列とreco列を分離
            std::vector<const ColumnSchema*> pred_cols, reco_cols;
            for (const auto& col : columns) {
                if (std::strncmp(col->name, "pred.", 5) == 0) {
                    pred_cols.push_back(col);
                } else if (std::strncmp(col->name, "reco.", 5) == 0) {
                    reco_cols.push_back(col);
                }
            }

            // 各Recordについて
            for (const auto& [_, r] : map) {
                // 1. pred部分を書き込み
                for (const auto& col : pred_cols) {
                    write_binary_value(out, col, r, nullptr);
                }

                // 2. reco数を書き込み
                uint32_t num_recos = static_cast<uint32_t>(r.recos.size());
                out.write(reinterpret_cast<const char*>(&num_recos), sizeof(num_recos));

                // 3. 各recoを書き込み
                for (const auto& reco : r.recos) {
                    for (const auto& col : reco_cols) {
                        write_binary_value(out, col, r, &reco);
                    }
                }
            }
        } else {
            // テキスト形式
            // まとめられた Record を出力
            for (const auto& [_, r] : map) {
                if (r.recos.empty()) {
                    write_line(out, r, {});
                } else {
                    for (const auto& rc : r.recos) {
                        write_line(out, r, &rc);
                    }
                }
            }
        }
    }

    void WriteAll(std::ostream& os, const std::vector<Record>& records,
               const std::vector<std::string>& custom_cols = {}) {
        if (!custom_cols.empty()) {
            SetColumns(custom_cols);
        }
        WriteHeader(os);
        WriteRecords(os, records);
    }

    void WriteAll(const std::string& filename, const std::vector<Record>& records,
               const std::vector<std::string>& custom_cols = {}) {
        if (!custom_cols.empty()) {
            SetColumns(custom_cols);
        }

        // 明示的な指定がない場合のみ、ファイル名からフォーマットを自動判定
        if (!format_explicit) SetFormatFromFilename(filename);

        if (format == "binary") {
            std::ofstream ofs(filename, std::ios::binary);
            if (!ofs) {
                throw std::runtime_error("Cannot open file: " + filename);
            }
            WriteHeader(ofs);
            WriteRecords(ofs, records);
        } else {
            std::ofstream ofs(filename);
            if (!ofs) {
                throw std::runtime_error("Cannot open file: " + filename);
            }
            WriteHeader(ofs);
            WriteRecords(ofs, records);
        }
    }

private:
    static const ColumnSchema* find_schema(const std::string& name) {
        for (const auto& s : SCHEMA) {
            if (name == s.name) return &s;
        }
        return {};
    }

    void write_line(std::ostream& out, const Record& r, const Reco* rc) const {
        for (size_t i = 0; i < columns.size(); ++i) {
            columns[i]->write(out, r, rc);
            if (i + 1 < columns.size()) out << delimiter;
        }
        out << "\n";
    }

    void write_binary_line(std::ostream& os, const Record& r, const Reco* rc) const {
        for (const auto& col : columns) {
            write_binary_value(os, col, r, rc);
        }
    }

    void write_binary_value(std::ostream& os, const ColumnSchema* col,
        const Record& r, const Reco* rc) const {
        switch (col->type) {
            case ColumnType::Int64: {
                int64_t val = col->get_int64(r, rc);
                os.write(reinterpret_cast<const char*>(&val), sizeof(val));
                break;
            }
            case ColumnType::Int: {
                int val = col->get_int(r, rc);
                os.write(reinterpret_cast<const char*>(&val), sizeof(val));
                break;
            }
            case ColumnType::Double: {
                double val = col->get_double(r, rc);
                os.write(reinterpret_cast<const char*>(&val), sizeof(val));
                break;
            }
            case ColumnType::Bool: {
                bool val = col->get_bool(r, rc);
                os.write(reinterpret_cast<const char*>(&val), sizeof(val));
                break;
            }
        }
    }

    // ファイル名からフォーマットを自動判定
    void SetFormatFromFilename(const std::string& filename) {
        if (filename.size() > 4) {
            std::string ext = filename.substr(filename.size() - 4);
            if (ext == ".tsv") SetFormat("tsv");
            else if (ext == ".bin") SetFormat("binary");
            else SetFormat("csv");
        }
    }

    void WriteBinary(std::ostream& os, const std::vector<Record>& records) const {
        // ヘッダー: マジックナンバーとバージョン
        os.write(reinterpret_cast<const char*>(&BINARY_MAGIC), sizeof(BINARY_MAGIC));
        os.write(reinterpret_cast<const char*>(&BINARY_VERSION), sizeof(BINARY_VERSION));

        // カラム情報
        uint32_t num_cols = static_cast<uint32_t>(columns.size());
        os.write(reinterpret_cast<const char*>(&num_cols), sizeof(num_cols));

        for (const auto& col : columns) {
            uint32_t name_len = static_cast<uint32_t>(std::strlen(col->name));
            os.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
            os.write(col->name, name_len);
            uint8_t type = static_cast<uint8_t>(col->type);
            os.write(reinterpret_cast<const char*>(&type), sizeof(type));
        }

        // データ: 各Recordのrecoごとに行を出力
        for (const auto& record : records) {
            if (record.recos.empty()) {
                write_binary_line(os, record, nullptr);
            } else {
                for (const auto& reco : record.recos) {
                    write_binary_line(os, record, &reco);
                }
            }
        }
    }
};

class FileReader {
public:
    std::string delimiter;
    std::string format;  // "csv", "tsv", "binary"
    bool format_explicit = false;  // 明示的にフォーマットが指定されたか

    // デフォルトコンストラクタ（CSVモード）
    FileReader() : delimiter(","), format("csv"), format_explicit(false) {}

    // コンストラクタ: フォーマットを指定（"csv", "tsv", または "binary"）
    FileReader(const std::string& fmt)
        : delimiter(","), format("csv"), format_explicit(true) {
        if (fmt == "tsv" || fmt == "\\t") {
            delimiter = "\t";
            format = "tsv";
        } else if (fmt == "binary" || fmt == "bin") {
            format = "binary";
        } else if (fmt == "csv" || fmt == ",") {
            // デフォルトのまま
        } else {
            throw std::runtime_error(
                "Invalid format specified: " + fmt);
        }
    }

    // フォーマットを設定
    void SetFormat(const std::string& fmt) {
        if (fmt == "tsv" || fmt == "\\t") {
            delimiter = "\t";
            format = "tsv";
        } else if (fmt == "binary" || fmt == "bin") {
            format = "binary";
        } else if (fmt == "csv" || fmt == ",") {
            delimiter = ",";
            format = "csv";
        } else {
            throw std::runtime_error(
                "Invalid format specified: " + fmt);
        }
        format_explicit = true;  // 明示的な指定をマーク
    }

    std::vector<Record> ReadAll(const std::string& filename) {
        // 明示的な指定がない場合のみ、ファイル名からフォーマットを自動判定
        if (!format_explicit) SetFormatFromFilename(filename);

        if (format == "binary") return ReadBinary(filename);
        else return ReadText(filename);
    }

private:
    std::vector<const ColumnSchema*> columns;

    void parse_columns_line(const std::string& s) {
        std::stringstream ss(s);
        std::string token;
        while (std::getline(ss, token, delimiter[0])) {
            auto pos = token.find('|');
            if (pos == std::string::npos) {
                throw std::runtime_error("Invalid column token: " + token);
            }
            std::string name = trim(token.substr(0, pos));
            const ColumnSchema* schema = find_schema(name);
            if (!schema) {
                throw std::runtime_error("Unknown column: " + name);
            }
            columns.push_back(schema);
        }
    }

    static std::string trim(std::string s) {
        auto is_space = [](char c) { return std::isspace(c); };
        s.erase(s.begin(), std::find_if_not(s.begin(), s.end(), is_space));
        s.erase(std::find_if_not(s.rbegin(), s.rend(), is_space).base(), s.end());
        return s;
    }

    static const ColumnSchema* find_schema(const std::string& name) {
        for (const auto& s : SCHEMA) {
            if (name == s.name) return &s;
        }
        return {};
    }

    // ファイル名からフォーマットを自動判定
    void SetFormatFromFilename(const std::string& filename) {
        if (filename.size() > 4) {
            std::string ext = filename.substr(filename.size() - 4);
            if (ext == ".tsv") SetFormat("tsv");
            else if (ext == ".bin") SetFormat("binary");
            else SetFormat("csv");
        }
    }

    std::vector<Record> ReadText(const std::string& filename) {
        // テキスト形式（CSV/TSV）
        std::ifstream ifs(filename);
        if (!ifs) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        // ヘッダーをパース
        std::vector<const ColumnSchema*> columns;
        bool has_chain_id = false;
        bool has_pred_pl = false;
        {
            std::string line;
            while (std::getline(ifs, line)) {
                if (line.rfind("# ", 0) == 0 && line.find("pred.") != std::string::npos) {
                    // カラム情報をパース
                    std::stringstream ss(line.substr(2));
                    std::string token;
                    while (std::getline(ss, token, delimiter[0])) {
                        auto pos = token.find('|');
                        if (pos == std::string::npos) {
                            throw std::runtime_error("Invalid column token: " + token);
                        }
                        std::string name = trim(token.substr(0, pos));
                        const ColumnSchema* schema = find_schema(name);
                        if (!schema) {
                            throw std::runtime_error("Unknown column: " + name);
                        }
                        columns.push_back(schema);
                        if (name == "pred.chain_id") has_chain_id = true;
                        if (name == "pred.pl") has_pred_pl = true;
                    }
                    break;
                }
            }
            if (columns.empty()) {
                throw std::runtime_error("Column header not found");
            }
            if (!has_chain_id || !has_pred_pl) {
                std::string missing = !has_chain_id && !has_pred_pl
                    ? "pred.chain_id, pred.pl"
                    : (!has_chain_id ? "pred.chain_id" : "pred.pl");
                throw std::runtime_error("Required columns missing: " + missing);
            }
        }

        // データ読み込み
        std::map<std::pair<int64_t, int>, Record> map;
        std::string line;
        RowData row;

        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::stringstream ss(line);
            std::string token;

            row = RowData{};
            Reco reco{};
            bool has_reco = false;

            for (size_t i = 0; i < columns.size(); ++i) {
                if (!std::getline(ss, token, delimiter[0])) break;
                token = trim(token);
                columns[i]->read(token, row, reco, has_reco);
            }

            if (has_reco) row.reco = reco;

            auto key = std::make_pair(row.pred.chain_id, row.pred.pl);
            auto& rec = map[key];
            rec.pred = row.pred;
            if (row.reco) {
                rec.recos.push_back(*row.reco);
            }
        }

        std::vector<Record> out;
        out.reserve(map.size());
        for (auto& [_, r] : map) {
            out.push_back(std::move(r));
        }
        return out;
    }

    std::vector<Record> ReadBinary(const std::string& filename) {
        // ファイルを開く
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        // マジックナンバーとバージョン確認
        uint32_t magic;
        uint16_t version;
        ifs.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        ifs.read(reinterpret_cast<char*>(&version), sizeof(version));

        if (magic != BINARY_MAGIC) {
            throw std::runtime_error("Invalid binary file format");
        }
        if (version != BINARY_VERSION) {
            throw std::runtime_error("Unsupported binary format version");
        }

        // カラム情報を読み込み
        uint32_t num_cols;
        ifs.read(reinterpret_cast<char*>(&num_cols), sizeof(num_cols));

        std::vector<const ColumnSchema*> columns;
        for (uint32_t i = 0; i < num_cols; ++i) {
            uint32_t name_len;
            ifs.read(reinterpret_cast<char*>(&name_len), sizeof(name_len));

            std::string col_name(name_len, '\0');
            ifs.read(&col_name[0], name_len);

            uint8_t type_byte;
            ifs.read(reinterpret_cast<char*>(&type_byte), sizeof(type_byte));

            const ColumnSchema* schema = find_schema(col_name);
            if (!schema) {
                throw std::runtime_error("Unknown column: " + col_name);
            }
            columns.push_back(schema);
        }

        // pred列とreco列を分離
        std::vector<const ColumnSchema*> pred_cols, reco_cols;
        for (const auto& col : columns) {
            if (std::strncmp(col->name, "pred.", 5) == 0) {
                pred_cols.push_back(col);
            } else if (std::strncmp(col->name, "reco.", 5) == 0) {
                reco_cols.push_back(col);
            }
        }

        // ストリーミング形式：EOFまで読み込み
        std::vector<Record> records;

        while (ifs.peek() != EOF) {
            if (!ifs) break;  // ストリーム エラーがあればループを抜ける

            Record record;
            RowData row{};
            Reco dummy_reco{};
            bool dummy_has_reco = false;

            // 1. pred部分を読み込み
            for (const auto& col : pred_cols) {
                switch (col->type) {
                    case ColumnType::Int64: {
                        int64_t val;
                        ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
                        if (!ifs) {
                            throw std::runtime_error(
                                "Failed to read Int64 column at record " + std::to_string(records.size())
                            );
                        }
                        if (col->set_int64) col->set_int64(val, row, dummy_reco, dummy_has_reco);
                        break;
                    }
                    case ColumnType::Int: {
                        int val;
                        ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
                        if (!ifs) {
                            throw std::runtime_error(
                                "Failed to read Int column at record " + std::to_string(records.size())
                            );
                        }
                        if (col->set_int) col->set_int(val, row, dummy_reco, dummy_has_reco);
                        break;
                    }
                    case ColumnType::Double: {
                        double val;
                        ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
                        if (!ifs) {
                            throw std::runtime_error(
                                "Failed to read Double column at record " + std::to_string(records.size())
                            );
                        }
                        if (col->set_double) col->set_double(val, row, dummy_reco, dummy_has_reco);
                        break;
                    }
                    case ColumnType::Bool: {
                        bool val;
                        ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
                        if (!ifs) {
                            throw std::runtime_error(
                                "Failed to read Bool column at record " + std::to_string(records.size())
                            );
                        }
                        if (col->set_bool) col->set_bool(val, row, dummy_reco);
                        break;
                    }
                }
            }
            record.pred = row.pred;

            // 2. reco数を読み込み
            uint32_t num_recos;
            ifs.read(reinterpret_cast<char*>(&num_recos), sizeof(num_recos));
            if (!ifs) {
                throw std::runtime_error("Failed to read reco count at record " + std::to_string(records.size()));
            }
            record.recos.reserve(num_recos);

            // 3. 各recoを読み込み
            for (uint32_t j = 0; j < num_recos; ++j) {
                Reco reco{};
                bool has_reco = true;

                for (const auto& col : reco_cols) {
                    switch (col->type) {
                        case ColumnType::Int64: {
                            int64_t val;
                            ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
                            if (!ifs) {
                                throw std::runtime_error(
                                    "Failed to read reco Int64 at record " + std::to_string(records.size())
                                    + ", reco " + std::to_string(j)
                                );
                            }
                            if (col->set_int64) col->set_int64(val, row, reco, has_reco);
                            break;
                        }
                        case ColumnType::Int: {
                            int val;
                            ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
                            if (!ifs) {
                                throw std::runtime_error(
                                    "Failed to read reco Int at record " + std::to_string(records.size())
                                    + ", reco " + std::to_string(j)
                                );
                            }
                            if (col->set_int) col->set_int(val, row, reco, has_reco);
                            break;
                        }
                        case ColumnType::Double: {
                            double val;
                            ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
                            if (!ifs) {
                                throw std::runtime_error(
                                    "Failed to read reco Double at record " + std::to_string(records.size())
                                    + ", reco " + std::to_string(j)
                                );
                            }
                            if (col->set_double) col->set_double(val, row, reco, has_reco);
                            break;
                        }
                        case ColumnType::Bool: {
                            bool val;
                            ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
                            if (!ifs) {
                                throw std::runtime_error(
                                    "Failed to read reco Bool at record " + std::to_string(records.size())
                                    + ", reco " + std::to_string(j)
                                );
                            }
                            if (col->set_bool) col->set_bool(val, row, reco);
                            break;
                        }
                    }
                }
                record.recos.push_back(reco);
            }

            records.push_back(std::move(record));
        }

        return records;
    }
};

} // namespace pred
