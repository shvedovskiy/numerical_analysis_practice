package ru.shvedovskiy.lib;

public class Matrix {
    public double[][] arr;
    public int n;
    public int m;
    public int rank;
    private int deletedColumns[];
    public int numberOfDeletedColumns;
    public double determinant;
    public static final double EPS = 1e-19;

    public Matrix(int m, int n) {
        this.n = n;
        this.m = m;
        arr = new double[m][n];
        deletedColumns = new int[n];
    }
    public Matrix(Matrix that) {
        n = that.n;
        m = that.m;
        arr = new double[m][n];
        deletedColumns = new int[n];
    }

    public Matrix add(Matrix that) {
        Matrix matrix = new Matrix(m, that.n);
        if (n == that.m) {
            for (int i = 0; i != m; ++i) {
                for (int j = 0; j < that.n; ++j) {
                    matrix.arr[i][j] = 0;
                    for (int k = 0; k != n; ++k) {
                        matrix.arr[i][j] += arr[i][k] * that.arr[k][j];
                    }
                }
            }
        }
        return matrix;
    }
    public Matrix subtract(Matrix that) {
        Matrix matrix = new Matrix(that);
        for (int i = 0; i != m; ++i) {
            for(int j = 0; j != n; ++j) {
                matrix.arr[i][j] = arr[i][j] - that.arr[i][j];
            }
        }
        return matrix;
    }
    public void extendedMatrix(Matrix first, Matrix second) {
        for (int i = 0; i != first.m; ++i) {
            for (int j = 0; j != first.n; ++j) {
                arr[i][j] = first.arr[i][j];
            }
        }
        for (int i = 0; i != second.m; ++i) {
            arr[i][first.n] = second.arr[i][0];
        }
        deletedColumns = new int[n];
    }

    private void deleteColumn(int column) {
        deletedColumns[column] = 1;
        ++numberOfDeletedColumns;

        Matrix tmp = new Matrix(m,n - 1);
        int numberOfDeleted;
        for (int i = 0; i != m; ++i) {
            numberOfDeleted = 0;
            for (int j = 0; j != n; ++j) {
                if (j != column) {
                    tmp.arr[i][j - numberOfDeleted] = arr[i][j];
                } else {
                    ++numberOfDeleted;
                }
            }
        }
        this.copyMatrix(tmp);
    }
    public void deleteZeroColumns() {
        int numberOfCols = n;
        int numberOfDeleted = 0;

        for (int j = 0; j != numberOfCols; ++j) {
            if (isNull(j - numberOfDeleted)) {
                deleteColumn(j - numberOfDeleted);
                ++numberOfDeleted;
            }
        }
    }
    public Matrix getColumn(int j) {
        Matrix res = new Matrix (m, 1);
        for (int i = 0; i != m; ++i)
        {
            res.arr[i][0] = arr[i][j];
        }
        return res;
    }
    public void setColumn(Matrix set_matrix, int j) {
        for (int i = 0; i != m; ++i) {
            arr[i][j] = set_matrix.arr[i][0];
        }
    }
    public double maxInColumn(int j) {
        double max = Math.abs(arr[0][j]);
        for (int i = 1; i != m; ++i) {
            if (max < Math.abs(arr[i][j])) {
                max = Math.abs(arr[i][j]);
            }
        }
        return max;
    }
    private void returnRow(int i) {
        Matrix H = new Matrix(this.m + 1, n);
        boolean isFound = false;

        for (int j = 0; j != this.m + 1; ++j) {
            if (i == j) {
                H.arr[j][0] = -5;
                isFound = true;
            } else {
                if (isFound) {
                    H.arr[j][0] = this.arr[j - 1][0];
                } else {
                    H.arr[j][0] = this.arr[j][0];
                }
            }
        }
        this.copyMatrix(H);
    }

    public double determinantOfTriangleMatrix() {
        determinant = 1;
        for (int i = 0; i != n; ++i) {
            if (Math.abs(arr[i][i]) < EPS) {
                return determinant = 0;
            } else {
                determinant *= arr[i][i];
            }
        }
        return determinant;
    }
    public int rankOfTriangleMatrix() {
        int res;

        if (n < m) {
            res = n;
        } else {
            res = m;
        }
        for (int i = 0; i != m; ++i) {
            if (arr[i][i] == 0) {
                res--;
            }
        }
        return res;
    }
    private int rankOfRectangleMatrix() {
        int res = m;
        int numberOfNulls;

        for (int i = 0; i != m; ++i) {
            numberOfNulls = 0;
            for (int j = 0; j != n; ++j) {
                if (arr[i][j] == 0) {
                    ++numberOfNulls;
                }
            }
            if (numberOfNulls >= n) {
                --res;
            }
        }
        return res;
    }
    public int rankOfMatrix() {
        int res;
        for (int k = 0; k != m; ++k) {
            for (int i = k + 1; i != m; ++i) {
                double factor = arr[i][k] / arr[k][k];

                if (arr[i][k] == 0 || Math.abs(arr[k][k]) < EPS) {
                    factor = 0;
                }
                for (int j = k; j != n; ++j) {
                    arr[i][j] -= factor * arr[k][j];
                }
            }
        }
        for (int i = 0; i != m; ++i) {
            for (int j = 0; j != n; ++j) {
                if (arr[i][j] <= EPS) {
                    arr[i][j] = Math.round(arr[i][j]);
                }
            }
        }
        res = this.rankOfRectangleMatrix();
        return res;
    }
    public double rate() {
        double res = 0;
        for (int i = 0; i != m; ++i) {
            for (int j = 0; j < n; ++j) {
                res += arr[i][j] * arr[i][j];
            }
        }
        res = Math.sqrt(res);
        return res;
    }

    public void setZeroMatrix() {
        for (int i = 0; i != m; ++i) {
            for (int j = 0; j != n; ++j) {
                arr[i][j] = 0;
            }
        }
    }
    public void createRandomMatrix() {
        for (int i = 0; i != m; ++i) {
            for (int j = 0; j != n; ++j) {
                arr[i][j] = (int) (Math.random() * 10);
            }
        }
    }
    public void getFullMatrix(Matrix that) {
        for (int i = 0; i != that.n; ++i) {
            if (that.deletedColumns[i] == 1) {
                this.returnRow(i);
            }
        }
    }
    public void setEMatrix() {
        for (int i = 0; i != m; ++i) {
            arr[i][i] = 1;
        }
    }
    private boolean isNull(int j) {
        for (int i = 0; i != m; ++i) {
            if (arr[i][j] != 0) {
                return false;
            }
        }
        return true;
    }

    public void copyMatrix(Matrix that) {
        n = that.n;
        m = that.m;
        arr = new double[m][n];

        for (int i = 0; i != m; ++i) {
            for (int j = 0; j != n; ++j) {
                arr[i][j] = that.arr[i][j];
            }
        }
        determinant = that.determinant;
        rank = that.rank;
    }
    public void printMatrix() {
        for (int i = 0; i != m; ++i) {
            for (int j = 0; j != n; ++j) {
                System.out.printf("%8.2f", arr[i][j]);
                System.out.print(" ");
            }
            System.out.println();
        }
        System.out.println();
    }
}
