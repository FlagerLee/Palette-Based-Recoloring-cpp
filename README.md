# Palette-Based-Recoloring-cpp

运行指令：

```[bash]
./main input_file_path num_pic num_color output_path
```

其中，input_file_path是给出的输入图片路径；num_pic是图片数量；num_color是颜色数量；output_path是
输出图片路径。示例：

```[bash]
./main video/%02d.png 21 5 test/
```

程序会生成一张Palette.png的调色板和Palette.txt的文件，里面储存RGB调色信息，可以直接用于输入。

```[bash]
./main input_file_path num_pic input_palette output_palette num_color output_path
```

其中，input_file_path是给出的输入图片路径，num_pic是图片数量。input_palette、output_palette是
存放palette的rgb信息的文件路径，其格式为每行一个rgb数值，用空格隔开。num_color是palette中
颜色数量，output_path是输出文件路径。示例：

```[bash]
./main video/%02d.png 21 test/Palette.txt test/OPalette.txt 5 test/
```
