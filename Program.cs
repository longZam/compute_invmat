using System.Numerics;
using System.Text;
using Spectre.Console;


static float ComputeDeterminant3x3(ReadOnlySpan<float> m)
{
    return m[0] * (m[4] * m[8] - m[5] * m[7])
                - m[1] * (m[3] * m[8] - m[5] * m[6])
                + m[2] * (m[3] * m[7] - m[4] * m[6]);
}


static Matrix4x4 ComputeCofactorMatrix(Matrix4x4 input)
{
    var result = new Matrix4x4();

    Span<float> submatrix = stackalloc float[9];

    int indexer = 0;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            // 소행렬 성분을 submatrix로
            for (int r = 0; r < 4; r++)
            {
                if (r == i) continue;

                for (int c = 0; c < 4; c++)
                {
                    if (c == j) continue;

                    submatrix[indexer++] = input[r, c];
                }
            }

            // i + j가 홀수면 부호 반전
            result[i, j] = (1 + ((i + j) % 2 * -2)) * ComputeDeterminant3x3(submatrix);

            indexer = 0;
        }
    }

    return result;
}


static bool TryComputeInverseByDeterminant(Matrix4x4 input, out Matrix4x4 output)
{
    output = Matrix4x4.Identity;

    var cofactors = ComputeCofactorMatrix(input);
    var determinant = 0f;

    for (int i = 0; i < 4; i++)
        determinant += input[0, i] * cofactors[0, i];

    if (determinant == 0)
        return false;

    output = Matrix4x4.Transpose(cofactors) * (1 / determinant);

    return true;
}


static bool TryComputeInverseByGaussJordan(Matrix4x4 input, out Matrix4x4 output)
{
    output = Matrix4x4.Identity;
    var inv = input;

    // 대각성분을 0이 아닌 값이 되도록 행 교환 기본행연산
    for (int i = 0; i < 4; i++)
    {
        if (Math.Abs(inv[i, i]) < float.Epsilon)
        {
            var big = i;

            for (int j = 0; j < 4; j++)
                if (Math.Abs(inv[j, i]) > Math.Abs(inv[big, i]))
                    big = j;
            
            if (big == i)
                return false;
            
            (inv[i], inv[big]) = (inv[big], inv[i]);
            (output[i], output[big]) = (output[big], output[i]);
        }
    }

    // 주대각선 아래 성분을 0으로 만드는 기본행연산 
    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            var multiplier = inv[j, i] / inv[i, i];
            inv[j] -= inv[i] * multiplier;
            output[j] -= output[i] * multiplier;
        }
    }

    // 주대각선 위 성분을 0으로 만드는 기본행연산
    for (int i = 3; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            var multiplier = inv[j, i] / inv[i, i];
            inv[j] -= inv[i] * multiplier;
            output[j] -= output[i] * multiplier;
        }
    }

    // 주대각선 성분을 1로 만드는 기본행연산
    for (int i = 0; i < 4; i++)
    {
        var divider = inv[i, i];
        inv[i] = inv[i] / divider;
        output[i] = output[i] / divider;
    }

    // 오차 보정
    for (int i = 0; i < 3; i++)
    {
        var error = Matrix4x4.Identity - input * output;
        output += output * error;
    }

    return true;
}


static string ToStringMatrix(int n, Matrix4x4 m)
{
    var stringBuilder = new StringBuilder();

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            stringBuilder.Append(m[i, j])
                         .Append(' ');
        }

        stringBuilder.AppendLine();
    }

    return stringBuilder.ToString();
}


static bool Equals(Matrix4x4 a, Matrix4x4 b, float epsilon)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            if (Math.Abs(a[i, j] - b[i, j]) > epsilon)
                return false;

    return true;
}


var n = await AnsiConsole.AskAsync<int>("정방행렬의 차수를 입력하세요: ");
var original = Matrix4x4.Identity;

// 4 미만의 차수를 가진 행렬을 블록 대각 행렬로 변환
for (int i = 0; i < n; i++)
{
    var line = await AnsiConsole.AskAsync<string>($"{i + 1}행: ");

    foreach ((var j, var value) in line.Split(' ').Select(float.Parse).Index())
        original[i, j] = value;
}


AnsiConsole.WriteLine();

AnsiConsole.WriteLine("행렬식으로 구한 역행렬:");
var sd = TryComputeInverseByDeterminant(original, out var byDeterminant);
AnsiConsole.WriteLine(sd ? ToStringMatrix(n, byDeterminant) : "역행렬이 존재하지 않습니다.");


AnsiConsole.WriteLine("가우스-조던 소거법으로 구한 역행렬:");
var sg = TryComputeInverseByGaussJordan(original, out var byGaussJordan);
AnsiConsole.WriteLine(sg ? ToStringMatrix(n, byGaussJordan) : "역행렬이 존재하지 않습니다.");

if (sd != sg)
{
    AnsiConsole.WriteLine("두 방법의 결과가 상이합니다.");
    return;
}

// 1e-4 오차 범위 
AnsiConsole.WriteLine(Equals(byDeterminant, byGaussJordan, 0.001f) ?
"두 방법의 결과가 동일합니다." :
"두 방법의 결과가 상이합니다.");