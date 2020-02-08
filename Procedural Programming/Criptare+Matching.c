#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#define presumed_detections 500

#if 1  // ACTIVATE / DEACTIVATE ALL CODE
typedef struct
{
	unsigned char b, g, r;
}pixel;

typedef struct
{
	unsigned int inaltime;
	unsigned int latime;
	pixel *vector_valori;
}BMP_image;	

typedef struct
{
	int indiceI;
	int indiceJ;
	double corelatie;
	pixel culoare;

}detectie;

typedef struct
{
	int inaltime;
	int latime;
	pixel **matrice;
}BMP_matrice;

//  |-----------------------------------|---/||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\
//  |   PARTEA INTAI - IMAGE CRIPTING   |---/||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\
//  |-----------------------------------|---/||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\
 
#if 1

BMP_image* citire_si_liniarizare(char *nume_img)
{
	BMP_image *img = (BMP_image*)malloc(sizeof(BMP_image) * 1);

	int i, j;
	FILE *fin = NULL;
	fin = fopen(nume_img, "rb");
	if (fin == NULL)
	{
		printf("nu am gasit imaginea sursa din care citesc\n");
		system("pause");
		exit(1);
	}

	fseek(fin, 18, SEEK_SET);
	fread(&img->latime, sizeof(unsigned int), 1, fin);
	fread(&img->inaltime, sizeof(unsigned int), 1, fin);
	printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n", img->latime, img->inaltime);

	char *temp_pixel = (char*)malloc(3 * sizeof(char));

	img->vector_valori = (pixel*)malloc((img->inaltime*img->latime) * sizeof(pixel));

	int padding;
	if (img->latime % 4 != 0)
		padding = 4 - (3 * img->latime) % 4;
	else
		padding = 0;

	fseek(fin, 54, SEEK_SET);

	for (i = img->inaltime - 1; i >= 0; i--)
	{
		for (j = 0; j < img->latime; j++)
		{
			fread(temp_pixel, sizeof(char), 3, fin);
			img->vector_valori[i*img->latime + j].b = temp_pixel[0];
			img->vector_valori[i*img->latime + j].g = temp_pixel[1];
			img->vector_valori[i*img->latime + j].r = temp_pixel[2];
			fflush(fin);
		}

		fseek(fin, padding, SEEK_CUR);
		
	}

	fclose(fin);
	free(temp_pixel);
	return img;
}

unsigned int* xorshift32(int seed, int size)
{
	unsigned int r, k;
	r = seed;
	unsigned int *output = (unsigned int*)malloc(size * sizeof(unsigned int));
	for (k = 0; k < size; k++)
	{
		r = r ^ r << 13;
		r = r ^ r >> 17;
		r = r ^ r << 5;
		output[k] = r;
		//pt ca populez de la 0 si nu de la 1,oriunde folosesc aleator, scad 1
	}
	return output;
}

BMP_image* permutare(BMP_image *img, unsigned int *aleatoare)
{
	pixel p;

	int *sigma = (int*)malloc((img->inaltime*img->latime) * sizeof(int));
	for (int i = 0; i <= img->inaltime*img->latime - 1; i++)
		sigma[i] = i;

	for (int i = img->latime*img->inaltime - 1, k = 0; i >= 1; i--)
	{
		unsigned int r = aleatoare[k++] % (i + 1);
		int aux = sigma[r];
		sigma[r] = sigma[i];
		sigma[i] = aux;
	}

	BMP_image *intermediare;
	intermediare = (BMP_image*)malloc(sizeof(BMP_image));

	intermediare->inaltime = img->inaltime;
	intermediare->latime = img->latime;

	intermediare->vector_valori = (pixel*)malloc((img->latime*img->inaltime) * sizeof(pixel));

	for (int i = 0; i <= img->inaltime*img->latime - 1; i++)
		intermediare->vector_valori[sigma[i]] = img->vector_valori[i];

	free(sigma);
	return intermediare;
}

void criptare(BMP_image *img, unsigned int *aleatoare, char *secret_key)
{
	FILE *fin = fopen(secret_key, "rt");
	if (fin == NULL)
	{
		fprintf(stderr, "nu s a gasit secret key");
		system("pause");
		exit(1);

	}
	unsigned int seed, SV;

	fscanf(fin, "%u", &seed);
	fscanf(fin, "%u", &SV);

	int i;
	img->vector_valori[0].b = SV ^ img->vector_valori[0].b ^ aleatoare[img->inaltime*img->latime - 1];
	img->vector_valori[0].g = (SV >> 8) ^ img->vector_valori[0].g ^ (aleatoare[img->inaltime*img->latime - 1] >> 8);
	img->vector_valori[0].r = (SV >> 16) ^ img->vector_valori[0].r ^ (aleatoare[img->inaltime*img->latime - 1] >> 16);

	for (i = 1; i <= img->inaltime*img->latime - 1; i++)
	{
		img->vector_valori[i].b = img->vector_valori[i - 1].b ^ img->vector_valori[i].b ^ aleatoare[img->inaltime*img->latime - 1 + i];//ma duc cu -1 mereu unde am aleatoriu
		img->vector_valori[i].g = img->vector_valori[i - 1].g ^ img->vector_valori[i].g ^ (aleatoare[img->inaltime*img->latime - 1 + i] >> 8);
		img->vector_valori[i].r = img->vector_valori[i - 1].r ^ img->vector_valori[i].r ^ (aleatoare[img->inaltime*img->latime - 1 + i] >> 16);
	}
	fclose(fin);
}

void creare_img(BMP_image *img, char *nume_img_sursa, char *nume_img_dest)
{
	FILE *fin = NULL;
	FILE *fout = NULL;
	fin = fopen(nume_img_sursa, "rb");
	fout = fopen(nume_img_dest, "wb");

	if (fin == NULL)
	{
		perror("fopen");
		fprintf(stderr, "can't open file %s", nume_img_sursa);
		exit(1);
	}

	if (fout == NULL)
	{
		perror("fopen");
		fprintf(stderr, "can't open file %s", nume_img_dest);
		exit(1);
	}

	unsigned char c;
	int i, j;

	int padding;
	if (img->latime % 4 != 0)
		padding = 4 - (3 * img->latime) % 4;
	else
		padding = 0;

	char header[54];
	fread(header, sizeof(char), 54, fin);
	fwrite(header, sizeof(char), 54, fout);

	for (i = img->inaltime - 1; i >= 0; i--)
	{
		for (j = 0; j < img->latime; j++)
		{
			fwrite(&img->vector_valori[i*img->latime + j].b, sizeof(char), 1, fout);
			fflush(fout);
			fwrite(&img->vector_valori[i*img->latime + j].g, sizeof(char), 1, fout);
			fflush(fout);
			fwrite(&img->vector_valori[i*img->latime + j].r, sizeof(char), 1, fout);
			fflush(fout);
		}
		if (i != 0)
			fseek(fout, padding, SEEK_CUR);
	}

	char *end_of_file = (char*)calloc(padding, sizeof(char));
	fwrite(end_of_file, sizeof(char), padding, fout);

	free(end_of_file);
	fclose(fin);
	fclose(fout);
}

void testul_chi_patrat(BMP_image *img)
{
	unsigned int *frecventa_blue = (int*)malloc(256 * sizeof(unsigned int));
	unsigned int *frecventa_green = (int*)malloc(256 * sizeof(unsigned int));
	unsigned int *frecventa_red = (int*)malloc(256 * sizeof(unsigned int));

	memset(frecventa_blue, 0, sizeof(int) * 256);
	memset(frecventa_green, 0, sizeof(int) * 256);
	memset(frecventa_red, 0, sizeof(int) * 256);

	int k, i;

	double suma_blue = 0, suma_green = 0, suma_red = 0;
	double frecventa_medie = (img->inaltime*img->latime) / 256.0;

	for (i = 0; i < img->latime*img->inaltime; i++)
	{
		frecventa_blue[(unsigned char)img->vector_valori[i].b]++;
		frecventa_green[(unsigned char)img->vector_valori[i].g]++;
		frecventa_red[(unsigned char)img->vector_valori[i].r]++;
	}

	for (i = 0; i < 256; i++)
	{
		suma_blue = suma_blue + ((frecventa_blue[i] - frecventa_medie)*(frecventa_blue[i] - frecventa_medie)) / frecventa_medie;
		suma_green = suma_green + ((frecventa_green[i] - frecventa_medie)*(frecventa_green[i] - frecventa_medie)) / frecventa_medie;
		suma_red = suma_red + ((frecventa_red[i] - frecventa_medie)*(frecventa_red[i] - frecventa_medie)) / frecventa_medie;
	}

	printf("Red = %.2f \nGreen = %.2f \nBlue = %.2f \n\n", suma_red, suma_green, suma_blue);

	free(frecventa_blue);
	free(frecventa_green);
	free(frecventa_red);

}

BMP_image* decriptare(BMP_image *img, unsigned int *aleatoare, char *secret_key)
{
	FILE *fin = fopen(secret_key, "rt");
	if (fin == NULL)
	{
		fprintf(stderr, "Nu s-a gasit fisierul\n");
		system("pause");
		exit(1);
	}
	unsigned int seed, SV;

	fscanf(fin, "%u", &seed);
	fscanf(fin, "%u", &SV);

	BMP_image *vector_decriptat = (BMP_image*)malloc(sizeof(BMP_image));
	vector_decriptat->vector_valori = (pixel*)malloc((img->inaltime * img->latime) * sizeof(pixel));

	vector_decriptat->inaltime = img->inaltime;
	vector_decriptat->latime = img->latime;

	vector_decriptat->vector_valori[0].b = SV ^ img->vector_valori[0].b ^ (aleatoare[img->inaltime * img->latime - 1]);
	vector_decriptat->vector_valori[0].g = (SV >> 8) ^ img->vector_valori[0].g ^ (aleatoare[img->inaltime * img->latime - 1] >> 8);
	vector_decriptat->vector_valori[0].r = (SV >> 16) ^ img->vector_valori[0].r ^ (aleatoare[img->inaltime * img->latime - 1] >> 16);

	for (int i = 1; i < img->inaltime * img->latime; i++)
	{
		vector_decriptat->vector_valori[i].b = img->vector_valori[i - 1].b ^ img->vector_valori[i].b ^ aleatoare[img->inaltime*img->latime - 1 + i];//ma duc cu -1 mereu unde am aleatoriu
		vector_decriptat->vector_valori[i].g = img->vector_valori[i - 1].g ^ img->vector_valori[i].g ^ (aleatoare[img->inaltime*img->latime - 1 + i] >> 8);
		vector_decriptat->vector_valori[i].r = img->vector_valori[i - 1].r ^ img->vector_valori[i].r ^ (aleatoare[img->inaltime*img->latime - 1 + i] >> 16);
	}

	unsigned int *sigma = (unsigned int*)malloc((img->inaltime*img->latime) * sizeof(unsigned int));

	for (int i = 0; i <= img->inaltime*img->latime - 1; i++)
		sigma[i] = i;

	for (int i = img->latime*img->inaltime - 1, k = 0; i >= 1; i--)
	{
		unsigned int r = aleatoare[k++] % (i + 1);
		int aux = sigma[r];
		sigma[r] = sigma[i];
		sigma[i] = aux;
	}
	BMP_image *img_decriptata = (BMP_image*)malloc(sizeof(BMP_image));
	img_decriptata->vector_valori = (pixel*)malloc((img->inaltime * img->latime) * sizeof(pixel));

	img_decriptata->inaltime = img->inaltime;
	img_decriptata->latime = img->latime;

	unsigned int *sigma_invers = (unsigned int*)malloc((img->inaltime*img->latime) * sizeof(unsigned int));
	for (long i = (long)(img->latime*img->inaltime) - 1; i >= 0; i--)
		sigma_invers[sigma[i]] = i;

	for (unsigned i = 0; i < img->latime*img->inaltime; i++)
	{
		img_decriptata->vector_valori[sigma_invers[i]] = vector_decriptat->vector_valori[i];
	}
	free(sigma);
	free(vector_decriptat->vector_valori);
	free(vector_decriptat);
	return img_decriptata;
}

#endif

//  |--------------------------------------|---/||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\
//  |   PARTEA A 2 A - TEMPLATE MATCHING   |---/||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\
//  |--------------------------------------|---/||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\

#if 1

void grayscale_image(char* nume_fisier_sursa, char* nume_fisier_destinatie)
{
	FILE *fin, *fout;
	unsigned int dim_img, latime_img, inaltime_img;
	unsigned char pRGB[3], header[54], aux;

	fin = fopen(nume_fisier_sursa, "rb");
	if (fin == NULL)
	{
		printf("Fisierul %s nu poate fi accesat/gasit\n", nume_fisier_sursa);
		system("pause");
		exit(1);
	}

	fout = fopen(nume_fisier_destinatie, "wb+");

	fseek(fin, 2, SEEK_SET);
	fread(&dim_img, sizeof(unsigned int), 1, fin);
	
	fseek(fin, 18, SEEK_SET);
	fread(&latime_img, sizeof(unsigned int), 1, fin);
	fread(&inaltime_img, sizeof(unsigned int), 1, fin);
	
	fseek(fin, 0, SEEK_SET);
	unsigned char c;
	while (fread(&c, 1, 1, fin) == 1)
	{
		fwrite(&c, 1, 1, fout);
		fflush(fout);
	}
	fclose(fin);

	int padding;
	if (latime_img % 4 != 0)
		padding = 4 - (3 * latime_img) % 4;
	else
		padding = 0;

	fseek(fout, 54, SEEK_SET);
	int i, j;
	for (i = 0; i < inaltime_img; i++)
	{
		for (j = 0; j < latime_img; j++)
		{
			fread(pRGB, 3, 1, fout);
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
			fseek(fout, -3, SEEK_CUR);
			fwrite(pRGB, 3, 1, fout);
			fflush(fout);
		}
		fseek(fout, padding, SEEK_CUR);
	}
	fclose(fout);
}

BMP_matrice* creare_matrice(char *nume_img_sursa)
{
	FILE *fin = NULL;
	FILE *fout = NULL;
	fin = fopen(nume_img_sursa, "rb");
	fout = fopen("test_citire.bmp", "wb");
	if (fin == NULL)
	{
		fprintf(stderr, "Imaginea %s nu a putut fi gasita", nume_img_sursa);
		system("pause");
		exit(1);
	}
	if (fout == NULL)
	{
		fprintf(stderr, "Imaginea %s nu a putut fi gasita", "test_citire.bmp");
		system("pause");
		exit(1);
	}

	BMP_matrice *img = (BMP_matrice*)malloc(sizeof(BMP_matrice));

	fseek(fin, 18, SEEK_SET);
	fread(&img->latime, sizeof(unsigned int), 1, fin);
	fread(&img->inaltime, sizeof(unsigned int), 1, fin);
	
	img->matrice = (pixel**)malloc(img->inaltime * sizeof(pixel*));
	for (int i = 0; i < img->inaltime; i++)
		img->matrice[i] = (pixel*)malloc(img->latime * sizeof(pixel));

	char *temp_pixel = (char*)malloc(3 * sizeof(char));

	fseek(fin, 0, SEEK_SET);
	char header[54];
	fread(header, sizeof(char), 54, fin);
	fwrite(header, sizeof(char), 54, fout);

	int padding;
	if (img->latime % 4 != 0)
		padding = 4 - (3 * img->latime) % 4;
	else
		padding = 0;

	fseek(fin, 54, SEEK_SET);
	unsigned char a = 0;
	for (int i = img->inaltime - 1; i >= 0; i--)
	{
		for (int j = 0; j < img->latime; j++)
		{
			fread(temp_pixel, sizeof(char), 3, fin);
			img->matrice[i][j].b = temp_pixel[0];
			img->matrice[i][j].g = temp_pixel[1];
			img->matrice[i][j].r = temp_pixel[2];
			fwrite(&img->matrice[i][j].b, sizeof(char), 1, fout);
			fwrite(&img->matrice[i][j].g, sizeof(char), 1, fout);
			fwrite(&img->matrice[i][j].r, sizeof(char), 1, fout);
			fflush(fout);
		}

		fseek(fin, padding, SEEK_CUR);
		fwrite(&a, sizeof(char), padding, fout);
	}

	fclose(fin);
	fclose(fout);
	free(temp_pixel);
	return img;
}

double deviatie_s(double S_mediu, BMP_matrice *sablon)
{
	double suma_s = 0;
	double n = (double)(sablon->latime * sablon->inaltime);

	for (int i = 0; i < sablon->inaltime; i++)
		for (int j = 0; j < sablon->latime; j++)
		{
			suma_s = suma_s + (sablon->matrice[i][j].r - S_mediu)*(sablon->matrice[i][j].r - S_mediu);
		}
	suma_s = sqrt(suma_s / (n - 1));
	return suma_s;
}

double deviatie_fi(int i, int j, double fi_mediu, BMP_matrice *sablon, BMP_matrice *img)
{
	double suma_fi = 0;
	double n = (double)(sablon->latime * sablon->inaltime);

	for (int k = i; k < i + sablon->inaltime; k++)
		for (int q = j; q < j + sablon->latime; q++)
			suma_fi = suma_fi + (img->matrice[k][q].r - fi_mediu)*(img->matrice[k][q].r - fi_mediu);
		
	suma_fi = sqrt(suma_fi / (n - 1));
	return suma_fi;
}

void contur_fereastra(BMP_matrice **img, int W, int H, detectie *detectii, int i)
{
	for (int k = detectii[i].indiceI; k < detectii[i].indiceI + H; k++)
	{
		(*img)->matrice[k][detectii[i].indiceJ].r = detectii[i].culoare.r;
		(*img)->matrice[k][detectii[i].indiceJ].g = detectii[i].culoare.g;
		(*img)->matrice[k][detectii[i].indiceJ].b = detectii[i].culoare.b;

		(*img)->matrice[k][detectii[i].indiceJ + W - 1].r = detectii[i].culoare.r;
		(*img)->matrice[k][detectii[i].indiceJ + W - 1].g = detectii[i].culoare.g;
		(*img)->matrice[k][detectii[i].indiceJ + W - 1].b = detectii[i].culoare.b;
	}
	
	for (int k = detectii[i].indiceJ; k < detectii[i].indiceJ + W; k++)
	{
		(*img)->matrice[detectii[i].indiceI][k].r = detectii[i].culoare.r;
		(*img)->matrice[detectii[i].indiceI][k].g = detectii[i].culoare.g;
		(*img)->matrice[detectii[i].indiceI][k].b = detectii[i].culoare.b;

		(*img)->matrice[detectii[i].indiceI + H - 1][k].r = detectii[i].culoare.r;
		(*img)->matrice[detectii[i].indiceI + H - 1][k].g = detectii[i].culoare.g;
		(*img)->matrice[detectii[i].indiceI + H - 1][k].b = detectii[i].culoare.b;
	}
}

double calculeaza_S_mediu(BMP_matrice *sablon, double *S_mediu)
{
	for (int i = 0; i < sablon->inaltime; i++)
		for (int j = 0; j < sablon->latime; j++)
			*S_mediu = *S_mediu + sablon->matrice[i][j].r;
	*S_mediu = *S_mediu / (double)(sablon->latime * sablon->inaltime);
	return *S_mediu;
}

double calculeaza_fi_mediu(BMP_matrice *img, BMP_matrice *sablon, double *fi_mediu, int i, int j)
{
	for (int k = i; k < i + sablon->inaltime; k++)
		for (int q = j; q < j + sablon->latime; q++)
			*fi_mediu = *fi_mediu + img->matrice[k][q].r;

	*fi_mediu = *fi_mediu / (sablon->latime*sablon->inaltime);
	return *fi_mediu;
}

double suprapunere(int DIindiceI, int DIindiceJ, int DJindiceI, int DJindiceJ, BMP_matrice *sablon)
{
	double arie_di, arie_dj, intersectie, reuniune;

	if (DIindiceI + sablon->inaltime <= DJindiceI || DIindiceJ + sablon->latime <= DJindiceJ)
		return 0;

	if (DJindiceI + sablon->inaltime <= DIindiceI || DIindiceJ + sablon->latime <= DIindiceJ)
		return 0;

	arie_di = sablon->inaltime * sablon->latime;
	arie_dj = arie_di;
	intersectie = (sablon->inaltime - abs(DIindiceI - DJindiceI))*(sablon->latime - abs(DIindiceJ - DJindiceJ));
	reuniune = arie_di + arie_dj - intersectie;
	return intersectie / reuniune;

}

void eliminare_non_maxim(int *size_detectii, BMP_matrice *sablon, detectie **detectii)
{
	int aux = 0, i, j;

	for (i = 0; i < (*size_detectii); i++)
	{
		if ((*detectii)[i].corelatie == 0)
			continue;
		aux++;
		for (j = i + 1; j < (*size_detectii); j++)
			if (suprapunere((*detectii)[i].indiceI, (*detectii)[i].indiceJ, (*detectii)[j].indiceI, (*detectii)[j].indiceJ, sablon) > 0.2)
				(*detectii)[i].corelatie = -1;
		
	}

	for (int i = 0; i < aux; i++)
		if ((*detectii)[i].corelatie == -1)
		{
			for (int k = i; k < aux - 1; k++)
				(*detectii)[k] = (*detectii)[k + 1];
			aux--;
			i--;
		}
	(*size_detectii) = aux;
}

int comparator(const void *c1, const void *c2)
{
	detectie *a = ((detectie*)c1);
	detectie* b = ((detectie*)c2);
	if (b->corelatie > a->corelatie)
		return 1;
	return -1;
}

void sorteaza(detectie *detectii, int size_detectii)
{
	qsort(detectii, size_detectii, sizeof(detectie), comparator);
}

BMP_matrice* matching_si_corelatie(char *nume_img_sursa_GS, char **vector_sabloane, float prag, detectie **detectii, int *size_detectii, char *nume_img_sursa)
{
	printf("Se proceseaza ");
	*size_detectii = 0;
	int test = 0;
	double dev_s;
	double dev_fi;
	double corelatie = 0;
	BMP_matrice *img = creare_matrice(nume_img_sursa_GS);
	BMP_matrice *sablon = creare_matrice(vector_sabloane[0]);
	int W = sablon->latime;
	int H = sablon->inaltime;

	double fi_mediu = 0;
	double S_mediu = 0;

	for (int nr_sablon = 0; nr_sablon < 10; nr_sablon++)
	{
		BMP_matrice *sablon = creare_matrice(vector_sabloane[nr_sablon]);

		printf(". ");
		S_mediu = calculeaza_S_mediu(sablon, &S_mediu);
		dev_s = deviatie_s(S_mediu, sablon);

		for (int i = 0; i < img->inaltime - sablon->inaltime; i++)
		{
			for (int j = 0; j < img->latime - sablon->latime; j++)
			{
				calculeaza_fi_mediu(img, sablon, &fi_mediu, i, j);
				dev_fi = deviatie_fi(i, j, fi_mediu, sablon, img);
				for (int k = i; k < i + sablon->inaltime; k++)
					for (int q = j; q < j + sablon->latime; q++)
						corelatie += (img->matrice[k][q].r - fi_mediu) * (sablon->matrice[k - i][q - j].r - S_mediu);

				corelatie /= (dev_fi*dev_s);
				corelatie = corelatie / (sablon->inaltime*sablon->latime);

				if (corelatie >= prag)
				{
					test++;

					if (*size_detectii == 0)
						(*detectii) = (detectie*)malloc(sizeof(detectie) * presumed_detections);
					else if ((*size_detectii) >= presumed_detections)
						(*detectii) = (detectie*)realloc((*detectii), sizeof(detectie)*(1 + *size_detectii));

					(*detectii)[*size_detectii].corelatie = corelatie;
					(*detectii)[*size_detectii].indiceI = i;
					(*detectii)[*size_detectii].indiceJ = j;

					switch (nr_sablon)
					{
					case 0:
					{
						(*detectii)[*size_detectii].culoare.r = 255;
						(*detectii)[*size_detectii].culoare.g = 0;
						(*detectii)[*size_detectii].culoare.b = 0;
					}
					break;

					case 1:
					{
						(*detectii)[*size_detectii].culoare.r = 255;
						(*detectii)[*size_detectii].culoare.g = 255;
						(*detectii)[*size_detectii].culoare.b = 0;
					}
					break;

					case 2:
					{
						(*detectii)[*size_detectii].culoare.r = 0;
						(*detectii)[*size_detectii].culoare.g = 255;
						(*detectii)[*size_detectii].culoare.b = 0;
					}
					break;

					case 3:
					{
						(*detectii)[*size_detectii].culoare.r = 0;
						(*detectii)[*size_detectii].culoare.g = 255;
						(*detectii)[*size_detectii].culoare.b = 255;
					}
					break;

					case 4:
					{
						(*detectii)[*size_detectii].culoare.r = 255;
						(*detectii)[*size_detectii].culoare.g = 0;
						(*detectii)[*size_detectii].culoare.b = 255;
					}
					break;

					case 5:
					{
						(*detectii)[*size_detectii].culoare.r = 0;
						(*detectii)[*size_detectii].culoare.g = 0;
						(*detectii)[*size_detectii].culoare.b = 255;
					}
					break;

					case 6:
					{
						(*detectii)[*size_detectii].culoare.r = 192;
						(*detectii)[*size_detectii].culoare.g = 192;
						(*detectii)[*size_detectii].culoare.b = 192;
					}
					break;

					case 7:
					{
						(*detectii)[*size_detectii].culoare.r = 255;
						(*detectii)[*size_detectii].culoare.g = 140;
						(*detectii)[*size_detectii].culoare.b = 0;
					}
					break;

					case 8:
					{
						(*detectii)[*size_detectii].culoare.r = 128;
						(*detectii)[*size_detectii].culoare.g = 0;
						(*detectii)[*size_detectii].culoare.b = 128;
					}
					break;

					case 9:
					{
						(*detectii)[*size_detectii].culoare.r = 128;
						(*detectii)[*size_detectii].culoare.g = 0;
						(*detectii)[*size_detectii].culoare.b = 0;
					}
					break;

					default:
						break;
					}
					(*size_detectii)++;
				}

				corelatie = 0;
			}
		}
	}
	BMP_matrice	*img_originala = creare_matrice(nume_img_sursa);
	
	sorteaza(*detectii, *size_detectii);

	eliminare_non_maxim(&(*size_detectii), sablon, (&(*detectii)));

	for (int i = 0; i < (*size_detectii); i++)
		contur_fereastra(&img_originala, W, H, (*detectii), i);

	free(sablon);
	free(*detectii);
	free(img);

	printf("\n");
	return img_originala;
}

void afisare_matrice(BMP_matrice *img, char *nume_img_sursa, char *nume_img_dest)
{
	FILE *fin = NULL;
	FILE *fout = NULL;
	fin = fopen(nume_img_sursa, "rb");
	fout = fopen(nume_img_dest, "wb");

	if (fin == NULL)
	{
		perror("fopen");
		fprintf(stderr, "can't open file %s", nume_img_sursa);
		exit(1);
	}

	if (fout == NULL)
	{
		perror("fopen");
		fprintf(stderr, "can't write in file %s", nume_img_dest);
		exit(1);
	}

	unsigned char c;
	int i, j;
	int padding;
	if (img->latime % 4 != 0)
		padding = 4 - (3 * img->latime) % 4;
	else
		padding = 0;

	char header[54];
	fread(header, sizeof(char), 54, fin);
	fwrite(header, sizeof(char), 54, fout);
	fclose(fin);
	
	for (i = img->inaltime - 1; i >= 0; i--)
	{
		for (j = 0; j < img->latime; j++)
		{
			fwrite(&img->matrice[i][j].b, sizeof(char), 1, fout);
			fflush(fout);
			fwrite(&img->matrice[i][j].g, sizeof(char), 1, fout);
			fflush(fout);
			fwrite(&img->matrice[i][j].r, sizeof(char), 1, fout);
			fflush(fout);
		}
		if (i != 0)
			fseek(fout, padding, SEEK_CUR);
	}

	char *end_of_file = (char*)calloc(padding, sizeof(char));
	fwrite(end_of_file, sizeof(char), padding, fout);

	free(end_of_file);

	fclose(fout);
}

#endif

#if 1
int main()
{
	printf("A intrat in main : \n\n");
	
	int MENIU;
	//   ---PARTEA INTAI---   //
	printf("Pe care parte a proiectului vrei sa o alegi? \n\n");
	printf("(1) CRIPTAREA UNEI IMAGINI\n");
	printf("(2) RECUNOASTEREA DE PATTERN-URI\n");
	printf("(3) PROGRAM COMPLET ( AMANDOUA )\n");
	scanf("%d", &MENIU);

	switch (MENIU)
	{
	case 1:
	{
		char *nume_img_sursa = (char*)malloc(100 * sizeof(char));
		char *nume_img_dest = (char*)malloc(100 * sizeof(char));
		char *nume_img_decriptata = (char*)malloc(100 * sizeof(char));
		char *secret_key = (char*)malloc(100 * sizeof(char));

		///    CITIRE SI LINIARIZARE   +   TEST CHI   ///

		printf("\nCum sa se numeasca imaginea sursa pe care vrei sa o criptezi : \n");
		scanf("%s", nume_img_sursa);

		BMP_image *imagine_sursa = citire_si_liniarizare(nume_img_sursa);
		printf("TESTUL CHI INAINTE DE CRIPTARE :\n\n");
		testul_chi_patrat(imagine_sursa);

		///   XORSHIFT32   ///

		printf("\nCare este numele fisierului care contine cheia secreta : \n");
		scanf("%s", secret_key);
		FILE *fin = fopen(secret_key, "rt");
		unsigned int seed;
		fscanf(fin, "%u", &seed);
		unsigned int *aleatoare = xorshift32(seed, imagine_sursa->latime*imagine_sursa->inaltime * 2 - 1);
		fclose(fin);

		///   PERRMUTARE   ///

		BMP_image *imagine_permutata = permutare(imagine_sursa, aleatoare);


		///   CRIPTARE    +   TEST CHI///

		printf("\nCum sa se numeasca imaginea criptata : \n");
		scanf("%s", nume_img_dest);

		criptare(imagine_permutata, aleatoare, secret_key);
		creare_img(imagine_permutata, nume_img_sursa, nume_img_dest);
		printf("TESTUL CHI DUPA CRIPTARE :\n\n");
		testul_chi_patrat(imagine_permutata);

		///   DECRIPTARE   ///

		printf("\nCum sa se numeasca imaginea decriptata : \n");
		scanf("%s", nume_img_decriptata);
		BMP_image *imagine_partial_decriptata = decriptare(imagine_permutata, aleatoare, secret_key);
		printf("\nSe creeaza imaginea decriptata...\n");
		creare_img(imagine_partial_decriptata, nume_img_sursa, nume_img_decriptata);
		printf("\nImaginea decriptata a fost creata\n");

		///   DEZALOCARE DE MEMORIE   ///

		free(nume_img_decriptata);
		free(nume_img_dest);
		free(nume_img_sursa);
		free(secret_key);
		free(imagine_permutata->vector_valori);
		free(imagine_permutata);
		free(imagine_sursa->vector_valori);
		free(imagine_sursa);
		free(aleatoare);

		break;
	}
	case 2:
	{
		char *imagine_sursa_1 = (char*)malloc(100 * sizeof(char));
		char *nume_img_finala = (char*)malloc(100 * sizeof(char));
		char *nume_sablon = (char*)malloc(100 * sizeof(char));
		char **vector_sabloane;
		int size_detectii = 0, nr_sabloane;
		float prag;
		BMP_matrice *img_matrice;
		detectie *detectii = NULL;

		vector_sabloane = (char**)malloc(10 * sizeof(char*));
		for (int i = 0; i < 10; i++)
			vector_sabloane[i] = (char*)malloc(100 * sizeof(char));

		printf("Cum se numeste imaginea sursa ? : ");
		scanf("%s", imagine_sursa_1);
		grayscale_image(imagine_sursa_1, "test_grayscaled.bmp");

		printf("Sablonul cifrei 0 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra0_gs.bmp");
		strcpy(vector_sabloane[0], "cifra0_gs.bmp");

		printf("Sablonul cifrei 1 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra1_gs.bmp");
		strcpy(vector_sabloane[1], "cifra1_gs.bmp");

		printf("Sablonul cifrei 2 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra2_gs.bmp");
		strcpy(vector_sabloane[2], "cifra2_gs.bmp");

		printf("Sablonul cifrei 3 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra3_gs.bmp");
		strcpy(vector_sabloane[3], "cifra3_gs.bmp");

		printf("Sablonul cifrei 4 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra4_gs.bmp");
		strcpy(vector_sabloane[4], "cifra4_gs.bmp");

		printf("Sablonul cifrei 5 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra5_gs.bmp");
		strcpy(vector_sabloane[5], "cifra5_gs.bmp");

		printf("Sablonul cifrei 6 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra6_gs.bmp");
		strcpy(vector_sabloane[6], "cifra6_gs.bmp");

		printf("Sablonul cifrei 7 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra7_gs.bmp");
		strcpy(vector_sabloane[7], "cifra7_gs.bmp");

		printf("Sablonul cifrei 8 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra8_gs.bmp");
		strcpy(vector_sabloane[8], "cifra8_gs.bmp");

		printf("Sablonul cifrei 9 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra9_gs.bmp");
		strcpy(vector_sabloane[9], "cifra9_gs.bmp");

		printf("Pragul dorit? : ");
		scanf("%f", &prag);

		img_matrice = matching_si_corelatie("test_grayscaled.bmp", vector_sabloane, prag, &detectii, &size_detectii, imagine_sursa_1);
		printf("Cum sa se numeasca imaginea finala? : ");
		scanf("%s", nume_img_finala);

		afisare_matrice(img_matrice, imagine_sursa_1, nume_img_finala);

		///   DEZALOCARE MEMORIE   ///

		free(imagine_sursa_1);
		free(nume_img_finala);
		free(nume_sablon);
		for (int i = 0; i < img_matrice->inaltime; i++)
			free(img_matrice->matrice[i]);
		free(img_matrice->matrice);
		free(img_matrice);

		for (int i = 0; i < 10; i++)
			free(vector_sabloane[i]);
		free(vector_sabloane);

		break;
	}
	case 3:
	{
		printf("PARTEA NR.1 INCEPE : \n\n");
		char *nume_img_sursa = (char*)malloc(100 * sizeof(char));
		char *nume_img_dest = (char*)malloc(100 * sizeof(char));
		char *nume_img_decriptata = (char*)malloc(100 * sizeof(char));
		char *secret_key = (char*)malloc(100 * sizeof(char));

		///    CITIRE SI LINIARIZARE   +   TEST CHI   ///

		printf("\nCum sa se numeasca imaginea sursa pe care vrei sa o criptezi : \n");
		scanf("%s", nume_img_sursa);

		BMP_image *imagine_sursa = citire_si_liniarizare(nume_img_sursa);
		printf("TESTUL CHI INAINTE DE CRIPTARE :\n\n");
		testul_chi_patrat(imagine_sursa);

		///   XORSHIFT32   ///

		printf("\nCare este numele fisierului care contine cheia secreta : \n");
		scanf("%s", secret_key);
		FILE *fin = fopen(secret_key, "rt");
		unsigned int seed;
		fscanf(fin, "%u", &seed);
		unsigned int *aleatoare = xorshift32(seed, imagine_sursa->latime*imagine_sursa->inaltime * 2 - 1);
		fclose(fin);

		///   PERRMUTARE   ///

		BMP_image *imagine_permutata = permutare(imagine_sursa, aleatoare);


		///   CRIPTARE    +   TEST CHI///

		printf("\nCum sa se numeasca imaginea criptata : \n");
		scanf("%s", nume_img_dest);

		criptare(imagine_permutata, aleatoare, secret_key);
		creare_img(imagine_permutata, nume_img_sursa, nume_img_dest);
		printf("TESTUL CHI DUPA CRIPTARE :\n\n");
		testul_chi_patrat(imagine_permutata);

		///   DECRIPTARE   ///

		printf("\nCum sa se numeasca imaginea decriptata : \n");
		scanf("%s", nume_img_decriptata);
		BMP_image *imagine_partial_decriptata = decriptare(imagine_permutata, aleatoare, secret_key);
		printf("\nSe creeaza imaginea decriptata...\n");
		creare_img(imagine_partial_decriptata, nume_img_sursa, nume_img_decriptata);
		printf("\nImaginea decriptata a fost creata\n");

		///   DEZALOCARE DE MEMORIE   ///

		free(nume_img_decriptata);
		free(nume_img_dest);
		free(nume_img_sursa);
		free(secret_key);
		free(imagine_permutata->vector_valori);
		free(imagine_permutata);
		free(imagine_sursa->vector_valori);
		free(imagine_sursa);
		free(aleatoare);

		printf("PARTEA NR.2 INCEPE : \n\n");

		char *imagine_sursa_1 = (char*)malloc(100 * sizeof(char));
		char *nume_img_finala = (char*)malloc(100 * sizeof(char));
		char *nume_sablon = (char*)malloc(100 * sizeof(char));
		char **vector_sabloane;
		int size_detectii = 0, nr_sabloane;
		float prag;
		BMP_matrice *img_matrice;
		detectie *detectii = NULL;

		vector_sabloane = (char**)malloc(10 * sizeof(char*));
		for (int i = 0; i < 10; i++)
			vector_sabloane[i] = (char*)malloc(100 * sizeof(char));

		printf("Cum se numeste imaginea sursa ? : ");
		scanf("%s", imagine_sursa_1);
		grayscale_image(imagine_sursa_1, "test_grayscaled.bmp");

		printf("Sablonul cifrei 0 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra0_gs.bmp");
		strcpy(vector_sabloane[0], "cifra0_gs.bmp");

		printf("Sablonul cifrei 1 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra1_gs.bmp");
		strcpy(vector_sabloane[1], "cifra1_gs.bmp");

		printf("Sablonul cifrei 2 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra2_gs.bmp");
		strcpy(vector_sabloane[2], "cifra2_gs.bmp");

		printf("Sablonul cifrei 3 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra3_gs.bmp");
		strcpy(vector_sabloane[3], "cifra3_gs.bmp");

		printf("Sablonul cifrei 4 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra4_gs.bmp");
		strcpy(vector_sabloane[4], "cifra4_gs.bmp");

		printf("Sablonul cifrei 5 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra5_gs.bmp");
		strcpy(vector_sabloane[5], "cifra5_gs.bmp");

		printf("Sablonul cifrei 6 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra6_gs.bmp");
		strcpy(vector_sabloane[6], "cifra6_gs.bmp");

		printf("Sablonul cifrei 7 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra7_gs.bmp");
		strcpy(vector_sabloane[7], "cifra7_gs.bmp");

		printf("Sablonul cifrei 8 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra8_gs.bmp");
		strcpy(vector_sabloane[8], "cifra8_gs.bmp");

		printf("Sablonul cifrei 9 se numeste : ");
		scanf("%s", nume_sablon);
		grayscale_image(nume_sablon, "cifra9_gs.bmp");
		strcpy(vector_sabloane[9], "cifra9_gs.bmp");

		printf("Pragul dorit? : ");
		scanf("%f", &prag);

		img_matrice = matching_si_corelatie("test_grayscaled.bmp", vector_sabloane, prag, &detectii, &size_detectii, imagine_sursa_1);
		printf("Cum sa se numeasca imaginea finala? : ");
		scanf("%s", nume_img_finala);

		afisare_matrice(img_matrice, imagine_sursa_1, nume_img_finala);

		///   DEZALOCARE MEMORIE   ///

		free(imagine_sursa_1);
		free(nume_img_finala);
		free(nume_sablon);
		for (int i = 0; i < img_matrice->inaltime; i++)
			free(img_matrice->matrice[i]);
		free(img_matrice->matrice);
		free(img_matrice);

		for (int i = 0; i < 10; i++)
			free(vector_sabloane[i]);
		free(vector_sabloane);

		break;
	}
	default:
	{
		printf("Aceasta optiune nu exista\n");
		break;
	}

	}
	
	system("pause");
}
#endif
#endif
