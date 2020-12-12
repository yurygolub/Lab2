using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab2
{
    class Vector3
    {
        private double[] arr = new double[3];

        public double[] Arr
        {
            get
            {
                return arr;
            }

            private set
            {
                arr = value;
            }
        }
        public Vector3()
        {
            for (int i = 0; i < arr.Length; i++) arr[i] = 0;
        }

        public Vector3(double x, double y, double z)
        {
            arr[0] = x;
            arr[1] = y;
            arr[2] = z;
        }

        public Vector3(Vector3 other)
        {
            arr[0] = other.arr[0];
            arr[1] = other.arr[1];
            arr[2] = other.arr[2];
        }

        public Vector3(Matrix other)
        {
            arr[0] = other.Arr[0, 0];
            arr[1] = other.Arr[1, 0];
            arr[2] = other.Arr[2, 0];
        }

        public static Vector3 Cross(Vector3 a1, Vector3 a2)
        {
            Vector3 result = new Vector3();
            result[0] = a1[1] * a2[2] - a1[2] * a2[1];
            result[1] = a1[2] * a2[0] - a1[0] * a2[2];
            result[2] = a1[0] * a2[1] - a1[1] * a2[0];
            return result;
        }

        public static double Dot(Vector3 a1, Vector3 a2)
        {
            double result = a1[0] * a2[0] + a1[1] * a2[1] + a1[2] * a2[2];
            return result;
        }

        public static double Norm(Vector3 a)
        {
            double result = 0;
            for (int i = 0; i < a.arr.Length; i++)
            {
                result += a[i] * a[i];
            }
            result = Math.Sqrt(result);
            return result;
        }

        public static Vector3 operator +(Vector3 vect, double value)
        {
            Vector3 temp = new Vector3();
            temp[0] = vect[0] + value;
            temp[1] = vect[1] + value;
            temp[2] = vect[2] + value;
            return temp;
        }

        public static Vector3 operator +(Vector3 vect1, Vector3 vect2)
        {
            Vector3 temp = new Vector3();
            temp[0] = vect1[0] + vect2[0];
            temp[1] = vect1[1] + vect2[1];
            temp[2] = vect1[2] + vect2[2];
            return temp;
        }

        public static Vector3 operator -(Vector3 vect, double value)
        {
            Vector3 temp = new Vector3();
            temp[0] = vect[0] - value;
            temp[1] = vect[1] - value;
            temp[2] = vect[2] - value;
            return temp;
        }

        public static Vector3 operator -(Vector3 vect1, Vector3 vect2)
        {
            Vector3 temp = new Vector3();
            temp[0] = vect1[0] - vect2[0];
            temp[1] = vect1[1] - vect2[1];
            temp[2] = vect1[2] - vect2[2];
            return temp;
        }

        public static Vector3 operator *(Vector3 vect, double value)
        {
            Vector3 temp = new Vector3();
            temp[0] = vect[0] * value;
            temp[1] = vect[1] * value;
            temp[2] = vect[2] * value;
            return temp;
        }

        public static Vector3 operator /(Vector3 vect, double value)
        {
            Vector3 temp = new Vector3();
            temp[0] = vect[0] / value;
            temp[1] = vect[1] / value;
            temp[2] = vect[2] / value;
            return temp;
        }

        public double this[int index]
        {
            get
            {
                return arr[index];
            }

            private set
            {
                arr[index] = value;
            }
        }

        public override string ToString()
        {
            return $"{arr[0]}\t{arr[1]}\t{arr[2]}";
        }
    }
}