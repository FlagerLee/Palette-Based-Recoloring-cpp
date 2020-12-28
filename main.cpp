#include "rbfcolor.h"
#include "png/lodepng.h"
#include <math.h>
#include <iostream>
#include <fstream>

using std::cin, std::cout, std::endl;
using std::ifstream;

int main(int argc, char** argv) {
    if(argc != 5 && argc != 7) return 0;
    if(argc == 5) {
        // file_path, num_pic, num_template, output_path
        string img_path(argv[1]);
        vector<unsigned char> pic;
        int pic_w, pic_h = 0;
        for(int i = 0; i < atoi(argv[2]); i ++) {
            char name[100];
            sprintf(name, img_path.c_str(), i);
            string img_file(name);
            vector<unsigned char> single_pic;
            single_pic.clear();
            
            unsigned int w, h;
            unsigned int error = lodepng::decode(single_pic, w, h, img_file.data());
            if(error != 0) {
                printf("Error when loading file %s\n", img_file.c_str());
                return 0;
            }
            pic_w = w; pic_h += h;
            pic.insert(pic.end(), single_pic.begin(), single_pic.end());
        }
        lodepng::encode("merge.png", pic, pic_w, pic_h);
        rbfcolor color(16, "merge.png");
        int num_palette = atoi(argv[3]);
        RGB* Palette = new RGB[8];
        color.gridacc_kmeans(num_palette, Palette);
        FILE* f = fopen((string(argv[4]) + "Palette.txt").c_str(), "w");
        for(int i = 0; i < num_palette; i ++) fprintf(f, "%.5lf %.5lf %.5lf\n", Palette[i].R, Palette[i].G, Palette[i].B);
        fclose(f);
        pic.clear();
        pic_w = 300, pic_h = 300 * num_palette;
        for(int i = 0; i < num_palette; i ++) {
            for(int j = 0; j < 90000; j ++) {
                unsigned char R = floor(Palette[i].R * 255.0);
                unsigned char G = floor(Palette[i].G * 255.0);
                unsigned char B = floor(Palette[i].B * 255.0);
                pic.push_back(R);
                pic.push_back(G);
                pic.push_back(B);
                pic.push_back(255);
            }
        }
        lodepng::encode(string(argv[4]) + "Palette.png", pic, pic_w, pic_h);
        delete []Palette;
    }
    else if(argc == 7) {
        // file_path, num_pic, input_palette, output_palette, num_palette, output_path
        int num_palette = atoi(argv[5]);
        RGB* input_palette = new RGB[8];
        RGB* output_palette = new RGB[8];
        RGB* delta_palette = new RGB[8];
        RGB* palette = new RGB[8];

        string input_palette_file(argv[3]);
        ifstream input_stream(input_palette_file);
        for(int i = 0; i < num_palette; i ++) {
            char buffer[100];
            input_stream.getline(buffer, 100);
            sscanf(buffer, "%lf %lf %lf", &input_palette[i].R, &input_palette[i].G, &input_palette[i].B);
        }
        input_stream.close();
        string output_palette_file(argv[4]);
        ifstream output_stream(output_palette_file);
        for(int i = 0; i < num_palette; i ++) {
            char buffer[100];
            output_stream.getline(buffer, 100);
            sscanf(buffer, "%lf %lf %lf", &output_palette[i].R, &output_palette[i].G, &output_palette[i].B);
        }

        int num_pic = atoi(argv[2]);
        for(int i = 0; i < num_pic; i ++) {
            delta_palette[i].R = (output_palette[i].R - input_palette[i].R) / (double)(num_pic - 1);
            delta_palette[i].G = (output_palette[i].G - input_palette[i].G) / (double)(num_pic - 1);
            delta_palette[i].B = (output_palette[i].B - input_palette[i].B) / (double)(num_pic - 1);
        }

        for(int i = 0; i < num_pic; i ++) {
            char name[100];
            sprintf(name, argv[1], i);
            string img_file(name);
            rbfcolor color(16, img_file);
            for(int j = 0; j < num_palette; j ++) {
                palette[j].R = input_palette[j].R + i * delta_palette[j].R;
                palette[j].G = input_palette[j].G + i * delta_palette[j].G;
                palette[j].B = input_palette[j].B + i * delta_palette[j].B;
            }
            color.draw(input_palette, palette, num_palette);
            sprintf(name, argv[6], i);
            string output_img(name);
            color.save_img(output_img);
        }
        delete []input_palette;
        delete []output_palette;
        delete []delta_palette;
        delete []palette;
        return 0;
    }
    return 0;
}

/*
int main(int argc, char** argv) {
    #ifndef script
    printf("Please type in image path: ");
    #endif
    string file_path;
    cin >> file_path;
    rbfcolor color(16, file_path);
    
    int num_palette = 5;
    RGB* Palette = new RGB[8];
    RGB* newPalette = new RGB[8];
    color.gridacc_kmeans(num_palette, Palette);
    for(int i = 0; i < num_palette; i ++) newPalette[i] = Palette[i];
    string input;
    while(std::getline(cin, input)) {
        if(input == "P" || input == "p") {
            // Palette
            #ifdef script
            for(int i = 0; i < num_palette; i ++) printf("%lf %lf %lf\n", newPalette[i].R, newPalette[i].G, newPalette[i].B);
            #else
            for(int i = 0; i < num_palette; i ++) printf("Palette %d: (RGB)(%d, %d, %d)\n", i, (int)floor(newPalette[i].R * 255.0), (int)floor(newPalette[i].G * 255.0), (int)floor(newPalette[i].B * 255.0));
            #endif
        }
        if(input == "C" || input == "c") {
            // Change Palette
            #ifndef script
            printf(">> Please type in Palette id: ");
            #endif
            int id; cin >> id;
            double input;
            #ifndef script
            printf(">> Please type in Palette R: ");
            #endif
            cin >> input;
            newPalette[id].R = input / 255.0;
            #ifndef script
            printf(">> Please type in Palette G: ");
            #endif
            cin >> input;
            newPalette[id].G = input / 255.0;
            #ifndef script
            printf(">> Please type in Palette B: ");
            #endif
            cin >> input;
            newPalette[id].B = input / 255.0;
        }
        if(input == "D" || input == "d") {
            // Draw
            #ifndef script
            printf(">> Please type in image path: ");
            #endif
            color.draw(Palette, newPalette, num_palette);
            string save_path;
            cin >> save_path;
            color.save_img(save_path);
        }
        if(input == "O" || input == "o") {
            // change origin
            #ifndef script
            printf(">> Please type in Palette id: ");
            #endif
            int id; cin >> id;
            double input;
            #ifndef script
            printf(">> Please type in Palette R: ");
            #endif
            cin >> input;
            Palette[id].R = input / 255.0;
            #ifndef script
            printf(">> Please type in Palette G: ");
            #endif
            cin >> input;
            Palette[id].G = input / 255.0;
            #ifndef script
            printf(">> Please type in Palette B: ");
            #endif
            cin >> input;
            Palette[id].B = input / 255.0;
        }
        if(input == "Q" || input == "q") break;
        if(input == "H" || input == "h") {
            // help
        }
        if(input == "V" || input == "v") {
            #ifndef script
            printf(">> Please type in save path: ");
            #endif
            string save_path;
            cin >> save_path;
            #ifndef script
            printf(">> Please type in number of pictures: ");
            #endif
            int num_pic;
            cin >> num_pic;
            RGB* delta = new RGB[8];
            RGB* delta_palette = new RGB[8];
            for(int i = 0; i < num_palette; i ++) {
                delta[i].R = (newPalette[i].R - Palette[i].R) / (double)(num_pic - 1);
                delta[i].G = (newPalette[i].G - Palette[i].G) / (double)(num_pic - 1);
                delta[i].B = (newPalette[i].B - Palette[i].B) / (double)(num_pic - 1);
            }
            for(int i = 0; i < num_pic; i ++) {
                for(int j = 0; j < num_palette; j ++) {
                    delta_palette[j].R = Palette[j].R + i * delta[j].R;
                    delta_palette[j].G = Palette[j].G + i * delta[j].G;
                    delta_palette[j].B = Palette[j].B + i * delta[j].B;
                }
                color.draw(Palette, delta_palette, num_palette);
                color.save_img(save_path + "output" + std::to_string(i) + ".png");
            }
            delete []delta;
            delete []delta_palette;
        }
        if(input == "R" || input == "r") {
            // reset
            for(int i = 0; i < num_palette; i ++) newPalette[i] = Palette[i];
        }
        if(input == "S" || input == "s") {
            // set num of palette
            #ifndef script
            printf(">> Please type in number of palette(2~8): ");
            #endif
            int new_num;
            cin >> new_num;
            if(new_num < 2 || new_num > 8) puts("Out Of Range!");
            else {
                num_palette = new_num;
                color.gridacc_kmeans(num_palette, Palette);
                for(int i = 0; i < num_palette; i ++) newPalette[i] = Palette[i];
            }
        }
    }
    return 0;
}
*/