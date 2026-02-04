/*
 * 文件名: modelselect.h
 * 文件作用: 模型选择对话框头文件
 * 功能描述:
 * 1. 定义模型选择对话框的类结构，包含UI指针和逻辑处理函数。
 * 2. 提供对外接口 getSelectedModelCode() 和 getSelectedModelName() 获取用户选择的模型结果。
 * 3. 提供 setCurrentModelCode() 接口支持根据模型ID反向初始化界面选项。
 * 修改记录:
 * 1. [调整] 根据最新的用户指令，重新定义了 Model 1-36 的分类映射逻辑。
 */

#ifndef MODELSELECT_H
#define MODELSELECT_H

#include <QDialog>

namespace Ui {
class ModelSelect;
}

class ModelSelect : public QDialog
{
    Q_OBJECT

public:
    /**
     * @brief 构造函数
     * @param parent 父窗口指针
     */
    explicit ModelSelect(QWidget *parent = nullptr);

    /**
     * @brief 析构函数
     */
    ~ModelSelect();

    /**
     * @brief 获取选中的模型代码
     * @return 返回模型ID字符串，例如 "modelwidget1"
     */
    QString getSelectedModelCode() const;

    /**
     * @brief 获取选中的模型显示名称
     * @return 返回模型名称字符串，例如 "夹层型储层试井解释1"
     */
    QString getSelectedModelName() const;

    /**
     * @brief 设置当前模型代码，用于打开窗口时回显之前的选择
     * @param code 模型ID字符串
     */
    void setCurrentModelCode(const QString& code);

private slots:
    /**
     * @brief 选项变更槽函数
     * 当任意下拉框选项发生变化时触发，重新计算对应的模型代码和名称
     */
    void onSelectionChanged();

    /**
     * @brief 确认按钮槽函数
     * 处理对话框确认逻辑
     */
    void onAccepted();

    /**
     * @brief 更新内外区选项槽函数
     * 根据当前选择的储层模型，动态更新内外区模型下拉框的内容
     */
    void updateInnerOuterOptions();

private:
    Ui::ModelSelect *ui;                // UI界面指针
    QString m_selectedModelCode;        // 当前选中的模型代码
    QString m_selectedModelName;        // 当前选中的模型名称

    /**
     * @brief 初始化下拉框选项
     * 设置井模型、储层模型、边界条件等下拉框的初始内容
     */
    void initOptions();

    bool m_isInitializing;              // 初始化标志位，防止信号递归触发
};

#endif // MODELSELECT_H
