using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab2
{
    class Program
    {
        private static double UTCToJulianDate(int year, int mounth, int day, int hour = 0, int min = 0, int sec = 0)
        {
            return 367 * year - (int)(7 * (year + (int)((mounth + 9) / 12.0)) / 4) + (int)(275 * mounth / 9.0) + 
                day + 1721013.5 + ((sec / 60.0 + min) / 60.0 + hour) / 24.0;
        }

        private static void Part1(double MU, double R, Vector3 r0, Vector3 V0)
        {
            Vector3 C = Vector3.Cross(r0, V0);
            Vector3 f = Vector3.Cross(V0, C) - r0 * MU / Vector3.Norm(r0);
            double h = Vector3.Norm(V0) * Vector3.Norm(V0) - 2 * MU / Vector3.Norm(r0);
            Vector3 ec = C / Vector3.Norm(C);
            double p = Vector3.Norm(C) * Vector3.Norm(C) / MU;
            double i = Math.Acos(ec[2]);
            double longitudeOfTheAscendingNode = Math.Asin(ec[0] / Math.Sin(i));
            double e = Math.Sqrt(1 + h * p / MU);
            double a = p / (1 - e * e);
            double trueAnomaly = Math.Asin(Vector3.Dot(V0, r0) / Vector3.Norm(r0) * Math.Sqrt(p / MU) / e);
            Vector3 tmp = new Vector3(Math.Cos(longitudeOfTheAscendingNode), Math.Sin(longitudeOfTheAscendingNode), 0);
            double latitudeArgument = Math.Acos(Vector3.Dot(tmp, r0) / Vector3.Norm(r0));
            double w0 = latitudeArgument - longitudeOfTheAscendingNode;
            double E0 = 2 * Math.Atan(Math.Sqrt((1 - e) / (1 + e)) * Math.Tan(trueAnomaly / 2));
            double Vmax = Math.Sqrt(MU / p) * (1 + e);
            double Vmin = Math.Sqrt(MU / p) * (1 - e);
            double T = 2 * Math.PI * Math.Sqrt(a * a * a / MU);
            double n = 2 * Math.PI / T;
            double rmin = p / (1 + e);
            double rmax = p / (1 - e);
            double Hmax = rmax - R;
            double Hmin = rmin - R;

            Console.WriteLine($"JulianDate:\t{UTCToJulianDate(2020, 3, 4, 16, 41, 0)}");

            Console.WriteLine($"a:\t{a}");
            Console.WriteLine($"p:\t{p}");
            Console.WriteLine($"i:\t{i}");
            Console.WriteLine($"e:\t{e}");
            Console.WriteLine($"longitudeOfTheAscendingNode:\t{longitudeOfTheAscendingNode}");
            Console.WriteLine($"w0:\t{w0}");
            Console.WriteLine($"trueAnomaly:\t{trueAnomaly}");
            Console.WriteLine($"Hmax:\t{Hmax}");
            Console.WriteLine($"Hmin:\t{Hmin}");
            Console.WriteLine($"T:\t{T}");
            Console.WriteLine($"n:\t{n}");
            Console.WriteLine($"Vmax:\t{Vmax}");
            Console.WriteLine($"Vmin:\t{Vmin}");
            Console.WriteLine($"latitudeArgument:\t{latitudeArgument}");

            Console.WriteLine($"h:\t{h}");
            Console.WriteLine($"C:\t{C}");
            Console.WriteLine($"f:\t{f}");
        }

        private static void Part2(double MU, double a, double omega, double e, double w, double i, double trueAnomaly)
        {
            e = e / 180 * Math.PI;
            omega = omega / 180 * Math.PI;
            w = w / 180 * Math.PI;
            i = i / 180 * Math.PI;
            trueAnomaly = trueAnomaly / 180 * Math.PI;
            double p = a * (1 - e * e);
            double sinTeta = Math.Sin(trueAnomaly);
            double cosTeta = Math.Cos(trueAnomaly);
            double ROrb = p / (1 + e * cosTeta);
            Vector3 rOrb = new Vector3(ROrb * cosTeta, ROrb * sinTeta, 0);
            Vector3 vOrb = new Vector3(-1 * sinTeta * Math.Sqrt(MU / p), (e + cosTeta) * Math.Sqrt(MU / p), 0);

            Vector3 r = new Vector3(Rotate3(-omega) * Rotate1(-i) * Rotate3(-w) * rOrb);
            Vector3 v = new Vector3(Rotate3(-omega) * Rotate1(-i) * Rotate3(-w) * vOrb);

            Console.WriteLine($"r:\n{r[0]}\n{r[1]}\n{r[2]}\n");
            Console.WriteLine($"v:\n{v[0]}\n{v[1]}\n{v[2]}\n");
        }

        private static void Part3(double MU, double R)
        {
            double eE = 0.081819221456;
            double Earth_angular_rate = 7.292115e-5;

            double epsilon = 1e-10;

            double phi_gd = (53 + 54 / 60.0 + 27 / 3600.0) / 180 * Math.PI;
            double longitude_gd = (27 + 33 / 60.0 + 52 / 3600.0) / 180 * Math.PI;
            double hellp = 0.23;

            double i = 97.4843;
            double omega = 290.6251;
            double e = 0.0015976;
            double w = 32.7870;
            double M0 = 41.0555;
            double n = 15.18422788105138;

            //double year = 2020;
            double day_tle = 265.12351617;

            double jd_epoch = UTCToJulianDate(2020, 1, 1, 0, 0, 0) + day_tle;
            double jd_begin = UTCToJulianDate(2020, 9, 21, 4, 0, 0);
            double jd_end = UTCToJulianDate(2020, 9, 21, 10, 0, 0);
            double delta_jd = 30.0 / 86400;
            int amount_points = (int)((jd_end - jd_begin) / delta_jd + 1);

            M0 = M0 / 180 * Math.PI;
            i = i / 180 * Math.PI;
            omega = omega / 180 * Math.PI;
            w = w / 180 * Math.PI;

            double a = Math.Pow(MU / Math.Pow(2 * Math.PI * n / 86400, 2), 1 / 3.0);
            double p = a * (1 - e * e);

            double Ce = R / Math.Sqrt(1 - eE * eE * Math.Sin(phi_gd) * Math.Sin(phi_gd));
            double Se = Ce * (1 - eE * eE);
            double r_delta = (Ce + hellp) * Math.Cos(phi_gd);
            double r_k = (Se + hellp) * Math.Sin(phi_gd);

            double Tut1 = (jd_begin - 2451545) / 36525;

            double GST1 = (100.4606184 + 36000.77005361 * Tut1 + 0.00038793 * Tut1 * Tut1 - 2.6e-8 * Math.Pow(Tut1, 3)) / 180 * Math.PI;

            double[] jd = new double[amount_points];
            double[] GST = new double[amount_points];
            double[] LST = new double[amount_points];
            double[] M = new double[amount_points];
            Vector3[] r_tso = new Vector3[amount_points];
            Vector3[] r_orb = new Vector3[amount_points];
            Vector3[] v_orb = new Vector3[amount_points];
            double[] r = new double[amount_points];
            double[] RA = new double[amount_points];
            double[] declination = new double[amount_points];
            Vector3[] r_eci = new Vector3[amount_points];
            Vector3[] v_eci = new Vector3[amount_points];
            Vector3[] v_ro_eci = new Vector3[amount_points];
            Vector3[] ro_eci = new Vector3[amount_points];
            Vector3[] ro_top = new Vector3[amount_points];
            Vector3[] v_top = new Vector3[amount_points];
            double[] ro_top_scs = new double[amount_points];
            double[] el = new double[amount_points];
            double[] beta = new double[amount_points];

            for (int index = 0; index < amount_points; index++)
            {
                jd[index] = jd_begin + (index - 1) * delta_jd;

                double UT1 = 86400 * (jd[index] - (int)jd[index]);

                GST[index] = GST1 + Earth_angular_rate * UT1;

                LST[index] = GST[index] + longitude_gd;

                r_tso[index] = new Vector3(r_delta * Math.Cos(LST[index]), r_delta * Math.Sin(LST[index]), r_k);

                M[index] = (M0 + n * (jd[index] - jd_epoch) * 2 * Math.PI) % (2 * Math.PI);

                double E_start;

                if (M[index] >= 0 && M[index] <= Math.PI)
                    E_start = M[index] + e;
                else
                    E_start = M[index] - e;

                int counter = 1;
                double E_previos = E_start;
                double E_next = E_previos + (M[index] + e * Math.Sin(E_previos) - E_previos) / (1 - e * Math.Cos(E_previos));

                while (Math.Abs(E_next - E_previos) > epsilon)
                {
                    E_previos = E_next;
                    E_next = E_previos + (M[index] - E_previos + e * Math.Sin(E_previos)) / (1 - e * Math.Cos(E_previos));
                    counter++;
                }

                double sinTeta = Math.Sqrt(1 - e * e) * Math.Sin(E_next) / (1 - e * Math.Cos(E_next));
                double cosTeta = (Math.Cos(E_next) - e) / (1 - e * Math.Cos(E_next));
                double R_orb = p / (1 + e * cosTeta);

                r_orb[index] = new Vector3(R_orb * cosTeta, R_orb * sinTeta, 0);
                v_orb[index] = new Vector3(-1 * sinTeta * Math.Sqrt(MU / p), (e + cosTeta) * Math.Sqrt(MU / p), 0);

                r[index] = Vector3.Norm(r_orb[index]);

                RA[index] = omega + Math.Atan2(Math.Cos(i) * sinTeta, cosTeta);

                if (RA[index] > 2 * Math.PI)
                    RA[index] = RA[index] - 2 * Math.PI;

                declination[index] = Math.Asin(Math.Sin(i) * sinTeta);

                r_eci[index] = new Vector3(Rotate3(-omega) * Rotate1(-i) * Rotate3(-w) * r_orb[index]);
                v_eci[index] = new Vector3(Rotate3(-omega) * Rotate1(-i) * Rotate3(-w) * v_orb[index]);
                Vector3 temp = new Vector3(0, 0, Earth_angular_rate);
                v_ro_eci[index] = new Vector3(v_eci[index] - Vector3.Cross(temp, r_eci[index]));
                ro_eci[index] = new Vector3(r_eci[index] - r_tso[index]);

                ro_top[index] = new Vector3(Rotate2(Math.PI / 2 - phi_gd) * Rotate3(LST[index]) * ro_eci[index]);
                v_top[index] = new Vector3(Rotate2(Math.PI / 2 - phi_gd) * Rotate3(LST[index]) * v_ro_eci[index]);

                ro_top_scs[index] = Vector3.Norm(ro_top[index]);
            }

            for (int index = 0; index < amount_points; index++)
            {
                RA[index] = RA[index] / Math.PI * 180;
                declination[index] = declination[index] / Math.PI * 180;

                el[index] = Math.Asin(ro_top[index].Arr[2] / ro_top_scs[index]) / Math.PI * 180;
                beta[index] = Math.Atan2(-ro_top[index].Arr[1], ro_top[index].Arr[0]) / Math.PI * 180;
            }

            PrintArray(jd, "jd");
            PrintArray(r_orb, "r_orb");
            PrintArray(v_orb, "v_orb");
            PrintArray(r_eci, "r_eci");
            PrintArray(v_eci, "v_eci");
            PrintArray(RA, "RA");
            PrintArray(declination, "declination");
            PrintArray(ro_top, "ro_top");
            PrintArray(v_top, "v_top");
            PrintArray(el, "el");
            PrintArray(beta, "beta");
            PrintArray(ro_top_scs, "ro_top_scs");

            SaveData(jd, "jd");
            SaveData(r_orb, "r_orb");
            SaveData(v_orb, "v_orb");
            SaveData(r_eci, "r_eci");
            SaveData(v_eci, "v_eci");
            SaveData(RA, "RA");
            SaveData(declination, "declination");
            SaveData(ro_top, "ro_top");
            SaveData(v_top, "v_top");
            SaveData(el, "el");
            SaveData(beta, "beta");
            SaveData(ro_top_scs, "ro_top_scs");
        }

        public static void SaveData<T>(T[] array, string title)
        {
            string path = $"{Environment.CurrentDirectory}\\{title}.txt";
            using (StreamWriter writer = File.CreateText(path))
            {
                string output = $"-----------------------{title}----------------------\n";
                
                for (int k = 0; k < array.Length; k++)
                {
                    output += array[k] + "\n";
                }
                writer.Write(output);
            }
        }

        public static void PrintArray<T>(T[] array, string title)
        {
            Console.WriteLine($"-----------------------{title}----------------------");
            for (int k = 0; k < array.Length; k++)
            {
                Console.WriteLine(array[k]);
            }
            Console.Write("\n\n\n\n\n\n");
        }

        private static Matrix Rotate1(double alpha)
        {
            double[,] temp = new double[3, 3]
            {
                {1, 0, 0 },
                {0, Math.Cos(alpha), Math.Sin(alpha)},
                {0, -Math.Sin(alpha),  Math.Cos(alpha)}
            };
            Matrix tempMatrix = new Matrix(temp);
            return tempMatrix;
        }

        private static Matrix Rotate2(double alpha)
        {
            double[,] temp = new double[3, 3]
            {
                {Math.Cos(alpha), 0, -Math.Sin(alpha) },
                {0, 1, 0},
                {Math.Sin(alpha), 0,  Math.Cos(alpha)}
            };
            Matrix tempMatrix = new Matrix(temp);
            return tempMatrix;
        }

        private static Matrix Rotate3(double alpha)
        {
            double[,] temp = new double[3, 3]
            {
                {Math.Cos(alpha), Math.Sin(alpha), 0 },
                {-Math.Sin(alpha), Math.Cos(alpha), 0 },
                {0, 0, 1 }
            };
            Matrix tempMatrix = new Matrix(temp);
            return tempMatrix;
        }

        static void Main(string[] args)
        {
            const double MU = 398600.4418, R = 6378.137;

            {
                Vector3 r0 = new Vector3(561.865, 4162.475, 5445.322);
                Vector3 V0 = new Vector3(0.825297, -6.057277, 4.5442540);
                Console.WriteLine("*********PART1*********\n");
                Part1(MU, R, r0, V0);
            }

            {
                const double a = 6884.1, omega = 92.2070, e = 0.0010320, w = 57.2909, i = 97.5014, trueAnomaly = 356.0200;
                Console.WriteLine("\n\n\n*********PART2*********");
                Part2(MU, a, omega, e, w, i, trueAnomaly);
            }

            {
                Console.WriteLine("\n\n\n*********PART3*********");
                Part3(MU, R);
            }
            Console.Read();
        }
    }
}