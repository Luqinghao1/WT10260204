/*
 * 文件名: modelsolver19_36.h
 * 文件作用: 压裂水平井页岩型复合模型核心计算类头文件
 * 功能描述:
 * 1. 提供 Model 19-36 (共18个) 模型的计算。
 * 2. 实现了Stehfest数值反演及压力导数计算。
 * 3. 核心算法采用瞬态平板模型(Transient Slab)来描述页岩型介质的流动特征。
 */

#ifndef MODELSOLVER19_36_H
#define MODELSOLVER19_36_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>
#include <QtConcurrent>

// 类型定义: <时间序列(t), 压力序列(Dp), 导数序列(Dp')>
using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver19_36
{
public:
    // 模型类型枚举 (对应 Model 19 - 36)
    enum ModelType {
        // --- 页岩型 + 页岩型 (19-24) ---
        Model_19 = 0, // 考虑井储表皮 + 无限大
        Model_20,     // 不考虑井储表皮 + 无限大
        Model_21,     // 考虑井储表皮 + 封闭
        Model_22,     // 不考虑井储表皮 + 封闭
        Model_23,     // 考虑井储表皮 + 定压
        Model_24,     // 不考虑井储表皮 + 定压

        // --- 页岩型 + 均质 (25-30) ---
        Model_25,     // 考虑井储表皮 + 无限大
        Model_26,     // 不考虑井储表皮 + 无限大
        Model_27,     // 考虑井储表皮 + 封闭
        Model_28,     // 不考虑井储表皮 + 封闭
        Model_29,     // 考虑井储表皮 + 定压
        Model_30,     // 不考虑井储表皮 + 定压

        // --- 页岩型 + 双重孔隙 (31-36) ---
        Model_31,     // 考虑井储表皮 + 无限大
        Model_32,     // 不考虑井储表皮 + 无限大
        Model_33,     // 考虑井储表皮 + 封闭
        Model_34,     // 不考虑井储表皮 + 封闭
        Model_35,     // 考虑井储表皮 + 定压
        Model_36      // 不考虑井储表皮 + 定压
    };

    explicit ModelSolver19_36(ModelType type);
    virtual ~ModelSolver19_36();

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

    // 边界元计算函数 (复用径向复合逻辑)
    double PWD_composite(double z, double fs1, double fs2, double M12, double LfD, double rmD, double reD,
                         int n_seg, int n_fracs, double spacingD, ModelType type);

    // 介质函数计算辅助
    double calc_fs_dual(double u, double omega, double lambda);       // 双重孔隙 f(s)

    // [新增] 页岩型 f(s) - 瞬态平板模型
    double calc_fs_shale(double u, double omega, double lambda);

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

#endif // MODELSOLVER19_36_H
