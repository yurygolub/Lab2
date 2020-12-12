using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab2
{
    class Matrix
    {
        private double[,] arr;

        public double[,] Arr
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

        public Matrix(int rows, int cols)
        {
            arr = new double[rows, cols];

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    arr[i, j] = 0;
        }

        public Matrix(double[,] _arr)
        {
            arr = new double[_arr.GetLength(0), _arr.GetLength(1)];
            for (int i = 0; i < arr.GetLength(0); i++)
                for (int j = 0; j < arr.GetLength(1); j++)
                    arr[i, j] = _arr[i, j];
        }

        public static Matrix operator *(Matrix matrix1, Matrix matrix2)
        {
            if (matrix1.arr.GetLength(1) != matrix2.arr.GetLength(0))
            {
                throw new Exception("Ошибка при умножении матриц");
            }

            Matrix temp = new Matrix(matrix1.arr.GetLength(0), matrix2.arr.GetLength(1));

            for (int i = 0; i < temp.arr.GetLength(0); i++)
            {
                for (int j = 0; j < temp.arr.GetLength(1); j++)
                {
                    for (int k = 0; k < matrix1.arr.GetLength(1); k++)
                    {
                        temp.arr[i, j] += matrix1.arr[i, k] * matrix2.arr[k, j];
                    }
                }
            }
            return temp;
        }

        public static Matrix operator *(Matrix matrix, Vector3 vector)
        {
            if (matrix.arr.GetLength(1) != vector.Arr.Length)
            {
                throw new Exception("Ошибка при умножении матриц");
            }

            Matrix temp = new Matrix(matrix.arr.GetLength(0), 1);

            for (int i = 0; i < temp.arr.GetLength(0); i++)
            {
                for (int k = 0; k < matrix.arr.GetLength(1); k++)
                {
                    temp.arr[i, 0] += matrix.arr[i, k] * vector.Arr[k];
                }
            }
            return temp;
        }

        public double this[int index1, int index2]
        {
            get
            {
                return arr[index1, index2];
            }

            private set
            {
                arr[index1, index2] = value;
            }
        }
    }
}