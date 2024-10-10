#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// 判断字符串是否为标准染色体号（完全匹配）
int is_valid_chromosome(const char *chrom) {
    if (strncmp(chrom, "chr", 3) != 0) {
        return 0; // 必须以 "chr" 开头
    }

    // 提取 "chr" 之后的部分
    const char *suffix = chrom + 3;

    // 判断是 1-22 之间的数字，且后面不能有额外字符
    char *endptr;
    long num = strtol(suffix, &endptr, 10); // 检查是否为数字
    if (endptr != suffix && *endptr == '\0' && num >= 1 && num <= 22) {
        return 1;
    }

    // 检查是否为 chrX, chrY 或 chrM
    if (strcmp(suffix, "X") == 0 || strcmp(suffix, "Y") == 0 || strcmp(suffix, "M") == 0) {
        return 1;
    }

    return 0; // 不符合标准染色体号
}

int process_input(FILE *input) {
    char line[1024]; // 假设每行最大长度为 1024 字节
    while (fgets(line, sizeof(line), input)) {
        // 因为 fgets 会保留换行符，所以我们需要处理换行符
        size_t len = strlen(line);
        if (len > 0 && line[len - 1] == '\n') {
            line[len - 1] = '\0'; // 去掉行尾的换行符
        }

        // 复制一份原始行数据，以防 strtok 修改原始行
        char line_copy[1024];
        strncpy(line_copy, line, sizeof(line_copy));
        
        // 用制表符分割第一列
        char *chrom = strtok(line, "\t"); 
        if (chrom != NULL && is_valid_chromosome(chrom)) {
            printf("%s\n", line_copy); // 输出原始的整行并加上换行符
        }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    // 如果有文件路径参数，则读取文件，否则从 stdin 读取（管道输入）
    if (argc == 2) {
        const char *file_path = argv[1];
        FILE *file = fopen(file_path, "r");
        if (file == NULL) {
            perror("Error opening file");
            return 1;
        }
        process_input(file);
        fclose(file);
    } else {
        process_input(stdin); // 从标准输入读取数据（适用于管道）
    }

    return 0;
}

// gcc -o chromosome_filter chromosome_filter.c