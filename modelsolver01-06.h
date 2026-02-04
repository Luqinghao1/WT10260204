/*
 * 文件名: modelsolver01-06.h
 * 文件作用: 压裂水平井复合模型(1-18)核心计算类头文件
 * 功能描述:
 * 1. 提供 Model 1-18 的试井模型计算。
 * 2. 核心计算逻辑基于双重孔隙和均质模型，根据用户定义映射为：
 * - 模型 1-6:  内区夹层型 + 外区夹层型 (计算逻辑：双孔+双孔)
 * - 模型 7-12: 径向复合均质 + 均质     (计算逻辑：均质+均质)
 * - 模型 13-18: 内区夹层型 + 外区均质   (计算逻辑：双孔+均质)
 * 3. 实现了Stehfest数值反演算法及压力导数计算。
 */

#ifndef MODELSOLVER01_06_H
#define MODELSOLVER01_06_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>
#include <QtConcurrent>

// 类型定义: <时间序列(t), 压力序列(Dp), 导数序列(Dp')>
using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver01_06
{
public:
    // 模型类型枚举
    enum ModelType {
        // --- 夹层型 + 夹层型 (Model 1-6) ---
        Model_1 = 0, // 考虑井储表皮 + 无限大
        Model_2,     // 不考虑井储表皮 + 无限大
        Model_3,     // 考虑井储表皮 + 封闭
        Model_4,     // 不考虑井储表皮 + 封闭
        Model_5,     // 考虑井储表皮 + 定压
        Model_6,     // 不考虑井储表皮 + 定压

        // --- 径向复合 (均质 + 均质) (Model 7-12) ---
        Model_7,     // 考虑井储表皮 + 无限大
        Model_8,     // 不考虑井储表皮 + 无限大
        Model_9,     // 考虑井储表皮 + 封闭
        Model_10,    // 不考虑井储表皮 + 封闭
        Model_11,    // 考虑井储表皮 + 定压
        Model_12,    // 不考虑井储表皮 + 定压

        // --- 夹层型 + 均质 (Model 13-18) ---
        Model_13,    // 考虑井储表皮 + 无限大
        Model_14,    // 不考虑井储表皮 + 无限大
        Model_15,    // 考虑井储表皮 + 封闭
        Model_16,    // 不考虑井储表皮 + 封闭
        Model_17,    // 考虑井储表皮 + 定压
        Model_18     // 不考虑井储表皮 + 定压
    };

    explicit ModelSolver01_06(ModelType type);
    virtual ~ModelSolver01_06();

    // 设置高精度计算模式
    void setHighPrecision(bool high);

    // 计算理论曲线接口
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime = QVector<double>());

    /**
     * @brief 获取模型名称
     * @param type 模型类型
     * @param verbose 是否显示详细信息
     */
    static QString getModelName(ModelType type, bool verbose = true);

    // 生成时间步
    static QVector<double> generateLogTimeSteps(int count, double startExp, double endExp);

private:
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                             std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                             QVector<double>& outPD, QVector<double>& outDeriv);

    // Laplace空间解主函数
    double flaplace_composite(double z, const QMap<QString, double>& p);

    // 边界元计算函数
    double PWD_composite(double z, double fs1, double fs2, double M12, double LfD, double rmD, double reD,
                         int n_seg, int n_fracs, double spacingD, ModelType type);

    // 数学辅助函数
    double scaled_besseli(int v, double x);
    double gauss15(std::function<double(double)> f, double a, double b);
    double adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth);

    // Stehfest算法辅助
    double getStehfestCoeff(int i, int N);
    void precomputeStehfestCoeffs(int N);
    double factorial(int n);

private:
    ModelType m_type;
    bool m_highPrecision;
    QVector<double> m_stehfestCoeffs;
    int m_currentN;
};

#endif // MODELSOLVER01_06_H
