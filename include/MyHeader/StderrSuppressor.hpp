#pragma once

#include <iostream>
#include <fstream>
#include <cstdio>
#include <io.h>

class StderrSuppressor {
public:
    void suppress() { // 標準エラー出力を抑制する
        // バックアップを保存
        old_cerr_buf = std::cerr.rdbuf();
        old_stderr_fd = _dup(_fileno(stderr));  // Cのstderrファイルディスクリプタを保存

        // C++標準エラー出力をNULにリダイレクト
        null_stream.open("NUL");
        std::cerr.rdbuf(null_stream.rdbuf());

        // CのstderrをNULにリダイレクト
        freopen("NUL", "w", stderr);
    }

    void restore() { // 標準エラー出力を復元する
        // C++出力の復元
        std::cerr.rdbuf(old_cerr_buf);

        // C出力の復元
        fflush(stderr);
        _dup2(old_stderr_fd, _fileno(stderr));
        _close(old_stderr_fd);

        null_stream.close();
    }

private:
    std::streambuf* old_cerr_buf = nullptr;
    int old_stderr_fd = -1;
    std::ofstream null_stream;
};