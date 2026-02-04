/*
 * 文件名: modelselect.cpp
 * 文件作用: 模型选择对话框逻辑实现
 * 功能描述:
 * 1. 初始化模型选择界面的各个下拉框选项。
 * 2. 实现储层模型与内外区模型选项的动态联动。
 * 3. 根据最新的模型对应表及用户指令，实现用户选择条件到模型ID (modelwidget1-36) 的映射。
 * - Model 1-6: 夹层型+夹层型
 * - Model 7-12: 径向复合 (均质+均质)
 * - Model 13-18: 夹层型+均质
 * - Model 19-36: 页岩型
 */

#include "modelselect.h"
#include "ui_modelselect.h"
#include <QDialogButtonBox>
#include <QPushButton>
#include <QDebug>

ModelSelect::ModelSelect(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ModelSelect),
    m_isInitializing(false)
{
    ui->setupUi(this);
    this->setStyleSheet("QWidget { color: black; font-family: Arial; }");

    initOptions();

    connect(ui->comboReservoirModel, SIGNAL(currentIndexChanged(int)), this, SLOT(updateInnerOuterOptions()));
    connect(ui->comboWellModel, SIGNAL(currentIndexChanged(int)), this, SLOT(onSelectionChanged()));
    connect(ui->comboReservoirModel, SIGNAL(currentIndexChanged(int)), this, SLOT(onSelectionChanged()));
    connect(ui->comboBoundary, SIGNAL(currentIndexChanged(int)), this, SLOT(onSelectionChanged()));
    connect(ui->comboStorage, SIGNAL(currentIndexChanged(int)), this, SLOT(onSelectionChanged()));
    connect(ui->comboInnerOuter, SIGNAL(currentIndexChanged(int)), this, SLOT(onSelectionChanged()));

    disconnect(ui->buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(ui->buttonBox, &QDialogButtonBox::accepted, this, &ModelSelect::onAccepted);
    connect(ui->buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    onSelectionChanged();
}

ModelSelect::~ModelSelect()
{
    delete ui;
}

void ModelSelect::initOptions()
{
    m_isInitializing = true;

    ui->comboWellModel->clear();
    ui->comboReservoirModel->clear();
    ui->comboBoundary->clear();
    ui->comboStorage->clear();
    ui->comboInnerOuter->clear();

    ui->comboWellModel->addItem("压裂水平井", "FracHorizontal");

    ui->comboReservoirModel->addItem("径向复合模型", "RadialComposite");         // 对应 Model 7-12 (均质+均质)
    ui->comboReservoirModel->addItem("夹层型径向复合模型", "InterlayerComposite"); // 对应 Model 1-6, 13-18
    ui->comboReservoirModel->addItem("页岩型径向复合模型", "ShaleComposite");      // 对应 Model 19-36
    ui->comboReservoirModel->addItem("混积型径向复合模型", "MixedComposite");      // 预留

    ui->comboBoundary->addItem("无限大外边界", "Infinite");
    ui->comboBoundary->addItem("封闭边界", "Closed");
    ui->comboBoundary->addItem("定压边界", "ConstantPressure");

    ui->comboStorage->addItem("考虑井储表皮", "Consider");
    ui->comboStorage->addItem("不考虑井储表皮", "Ignore");

    ui->comboWellModel->setCurrentIndex(0);
    ui->comboReservoirModel->setCurrentIndex(0);
    ui->comboBoundary->setCurrentIndex(0);
    ui->comboStorage->setCurrentIndex(0);

    m_isInitializing = false;
    updateInnerOuterOptions();
}

void ModelSelect::updateInnerOuterOptions()
{
    bool oldState = ui->comboInnerOuter->blockSignals(true);
    ui->comboInnerOuter->clear();
    QString currentRes = ui->comboReservoirModel->currentData().toString();

    if (currentRes == "RadialComposite") {
        // Model 7-12: 径向复合 -> 均质+均质
        ui->comboInnerOuter->addItem("均质+均质", "Homo_Homo");
    }
    else if (currentRes == "InterlayerComposite") {
        // Model 1-6: 夹层型+夹层型
        ui->comboInnerOuter->addItem("夹层型+夹层型", "Interlayer_Interlayer");
        // Model 13-18: 夹层型+均质
        ui->comboInnerOuter->addItem("夹层型+均质", "Interlayer_Homo");
    }
    else if (currentRes == "ShaleComposite") {
        // Model 19-36
        ui->comboInnerOuter->addItem("页岩型+页岩型", "Shale_Shale");
        ui->comboInnerOuter->addItem("页岩型+均质", "Shale_Homo");
        ui->comboInnerOuter->addItem("页岩型+双重孔隙", "Shale_Dual");
    }
    else if (currentRes == "MixedComposite") {
        ui->comboInnerOuter->addItem("混积型+混积型", "Mixed_Mixed");
        ui->comboInnerOuter->addItem("混积型+均质", "Mixed_Homo");
        ui->comboInnerOuter->addItem("混积型+双重孔隙", "Mixed_Dual");
    }

    if (ui->comboInnerOuter->count() > 0) {
        ui->comboInnerOuter->setCurrentIndex(0);
    }

    ui->comboInnerOuter->blockSignals(oldState);
    ui->label_InnerOuter->setVisible(true);
    ui->comboInnerOuter->setVisible(true);
}

void ModelSelect::setCurrentModelCode(const QString& code)
{
    m_isInitializing = true;
    QString numStr = code;
    numStr.remove("modelwidget");
    int id = numStr.toInt();

    if (id >= 1) {
        int idxWell = ui->comboWellModel->findData("FracHorizontal");
        if (idxWell >= 0) ui->comboWellModel->setCurrentIndex(idxWell);

        QString resData, ioData;

        if (id >= 1 && id <= 6) {
            // 夹层型 (1-6) -> 夹层+夹层
            resData = "InterlayerComposite";
            ioData = "Interlayer_Interlayer";
        }
        else if (id >= 7 && id <= 12) {
            // 径向复合 (7-12) -> 均质+均质
            resData = "RadialComposite";
            ioData = "Homo_Homo";
        }
        else if (id >= 13 && id <= 18) {
            // 夹层型 (13-18) -> 夹层+均质
            resData = "InterlayerComposite";
            ioData = "Interlayer_Homo";
        }
        else if (id >= 19 && id <= 36) {
            resData = "ShaleComposite";
            if (id <= 24) ioData = "Shale_Shale";
            else if (id <= 30) ioData = "Shale_Homo";
            else ioData = "Shale_Dual";
        }

        int idxRes = ui->comboReservoirModel->findData(resData);
        if (idxRes >= 0) {
            ui->comboReservoirModel->setCurrentIndex(idxRes);
            updateInnerOuterOptions();
        }

        QString bndData;
        int rem = (id - 1) % 6;
        if (rem == 0 || rem == 1) bndData = "Infinite";
        else if (rem == 2 || rem == 3) bndData = "Closed";
        else bndData = "ConstantPressure";

        int idxBnd = ui->comboBoundary->findData(bndData);
        if (idxBnd >= 0) ui->comboBoundary->setCurrentIndex(idxBnd);

        QString storeData = (id % 2 != 0) ? "Consider" : "Ignore";
        int idxStore = ui->comboStorage->findData(storeData);
        if (idxStore >= 0) ui->comboStorage->setCurrentIndex(idxStore);

        int idxIo = ui->comboInnerOuter->findData(ioData);
        if (idxIo >= 0) ui->comboInnerOuter->setCurrentIndex(idxIo);
    }

    m_isInitializing = false;
    onSelectionChanged();
}

void ModelSelect::onSelectionChanged()
{
    if (m_isInitializing) return;

    QString well = ui->comboWellModel->currentData().toString();
    QString res = ui->comboReservoirModel->currentData().toString();
    QString bnd = ui->comboBoundary->currentData().toString();
    QString store = ui->comboStorage->currentData().toString();
    QString io = ui->comboInnerOuter->currentData().toString();

    m_selectedModelCode = "";
    m_selectedModelName = "";

    // === 1. 径向复合模型 (Model 7-12) ===
    if (well == "FracHorizontal" && res == "RadialComposite") {
        if (io == "Homo_Homo") { // 7-12
            if (bnd == "Infinite") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget7"; m_selectedModelName = "径向复合储层试井解释7"; }
                else                     { m_selectedModelCode = "modelwidget8"; m_selectedModelName = "径向复合储层试井解释8"; }
            } else if (bnd == "Closed") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget9"; m_selectedModelName = "径向复合储层试井解释9"; }
                else                     { m_selectedModelCode = "modelwidget10"; m_selectedModelName = "径向复合储层试井解释10"; }
            } else if (bnd == "ConstantPressure") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget11"; m_selectedModelName = "径向复合储层试井解释11"; }
                else                     { m_selectedModelCode = "modelwidget12"; m_selectedModelName = "径向复合储层试井解释12"; }
            }
        }
    }
    // === 2. 夹层型径向复合模型 (Model 1-6, 13-18) ===
    else if (well == "FracHorizontal" && res == "InterlayerComposite") {
        if (io == "Interlayer_Interlayer") { // 1-6
            if (bnd == "Infinite") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget1"; m_selectedModelName = "夹层型储层试井解释1"; }
                else                     { m_selectedModelCode = "modelwidget2"; m_selectedModelName = "夹层型储层试井解释2"; }
            } else if (bnd == "Closed") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget3"; m_selectedModelName = "夹层型储层试井解释3"; }
                else                     { m_selectedModelCode = "modelwidget4"; m_selectedModelName = "夹层型储层试井解释4"; }
            } else if (bnd == "ConstantPressure") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget5"; m_selectedModelName = "夹层型储层试井解释5"; }
                else                     { m_selectedModelCode = "modelwidget6"; m_selectedModelName = "夹层型储层试井解释6"; }
            }
        }
        else if (io == "Interlayer_Homo") { // 13-18 (夹层型+均质)
            if (bnd == "Infinite") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget13"; m_selectedModelName = "夹层型储层试井解释13"; }
                else                     { m_selectedModelCode = "modelwidget14"; m_selectedModelName = "夹层型储层试井解释14"; }
            } else if (bnd == "Closed") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget15"; m_selectedModelName = "夹层型储层试井解释15"; }
                else                     { m_selectedModelCode = "modelwidget16"; m_selectedModelName = "夹层型储层试井解释16"; }
            } else if (bnd == "ConstantPressure") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget17"; m_selectedModelName = "夹层型储层试井解释17"; }
                else                     { m_selectedModelCode = "modelwidget18"; m_selectedModelName = "夹层型储层试井解释18"; }
            }
        }
    }
    // === 3. 页岩型径向复合模型 (Model 19-36) ===
    else if (well == "FracHorizontal" && res == "ShaleComposite") {
        if (io == "Shale_Shale") { // 19-24
            if (bnd == "Infinite") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget19"; m_selectedModelName = "页岩型储层试井解释1"; }
                else                     { m_selectedModelCode = "modelwidget20"; m_selectedModelName = "页岩型储层试井解释2"; }
            } else if (bnd == "Closed") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget21"; m_selectedModelName = "页岩型储层试井解释3"; }
                else                     { m_selectedModelCode = "modelwidget22"; m_selectedModelName = "页岩型储层试井解释4"; }
            } else if (bnd == "ConstantPressure") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget23"; m_selectedModelName = "页岩型储层试井解释5"; }
                else                     { m_selectedModelCode = "modelwidget24"; m_selectedModelName = "页岩型储层试井解释6"; }
            }
        }
        else if (io == "Shale_Homo") { // 25-30
            if (bnd == "Infinite") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget25"; m_selectedModelName = "页岩型储层试井解释7"; }
                else                     { m_selectedModelCode = "modelwidget26"; m_selectedModelName = "页岩型储层试井解释8"; }
            } else if (bnd == "Closed") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget27"; m_selectedModelName = "页岩型储层试井解释9"; }
                else                     { m_selectedModelCode = "modelwidget28"; m_selectedModelName = "页岩型储层试井解释10"; }
            } else if (bnd == "ConstantPressure") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget29"; m_selectedModelName = "页岩型储层试井解释11"; }
                else                     { m_selectedModelCode = "modelwidget30"; m_selectedModelName = "页岩型储层试井解释12"; }
            }
        }
        else if (io == "Shale_Dual") { // 31-36
            if (bnd == "Infinite") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget31"; m_selectedModelName = "页岩型储层试井解释13"; }
                else                     { m_selectedModelCode = "modelwidget32"; m_selectedModelName = "页岩型储层试井解释14"; }
            } else if (bnd == "Closed") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget33"; m_selectedModelName = "页岩型储层试井解释15"; }
                else                     { m_selectedModelCode = "modelwidget34"; m_selectedModelName = "页岩型储层试井解释16"; }
            } else if (bnd == "ConstantPressure") {
                if (store == "Consider") { m_selectedModelCode = "modelwidget35"; m_selectedModelName = "页岩型储层试井解释17"; }
                else                     { m_selectedModelCode = "modelwidget36"; m_selectedModelName = "页岩型储层试井解释18"; }
            }
        }
    }

    bool isValid = !m_selectedModelCode.isEmpty();

    if (isValid) {
        ui->label_ModelName->setText(m_selectedModelName);
        ui->label_ModelName->setStyleSheet("color: black; font-weight: bold; font-size: 14px;");
    } else {
        ui->label_ModelName->setText("该组合暂无已实现模型");
        ui->label_ModelName->setStyleSheet("color: red; font-weight: normal;");
    }

    QPushButton* okBtn = ui->buttonBox->button(QDialogButtonBox::Ok);
    if(okBtn) okBtn->setEnabled(isValid);
}

void ModelSelect::onAccepted() {
    if (!m_selectedModelCode.isEmpty()) accept();
}

QString ModelSelect::getSelectedModelCode() const { return m_selectedModelCode; }
QString ModelSelect::getSelectedModelName() const { return m_selectedModelName; }
