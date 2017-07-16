namespace SacBunt
{
    partial class MainFrm
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.textBoxResult = new System.Windows.Forms.TextBox();
            this.simulateBtn = new System.Windows.Forms.Button();
            this.SuspendLayout();
            // 
            // textBoxResult
            // 
            this.textBoxResult.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.textBoxResult.Location = new System.Drawing.Point(13, 41);
            this.textBoxResult.Multiline = true;
            this.textBoxResult.Name = "textBoxResult";
            this.textBoxResult.ScrollBars = System.Windows.Forms.ScrollBars.Vertical;
            this.textBoxResult.Size = new System.Drawing.Size(555, 382);
            this.textBoxResult.TabIndex = 12;
            // 
            // simulateBtn
            // 
            this.simulateBtn.Location = new System.Drawing.Point(12, 12);
            this.simulateBtn.Name = "simulateBtn";
            this.simulateBtn.Size = new System.Drawing.Size(553, 23);
            this.simulateBtn.TabIndex = 14;
            this.simulateBtn.Text = "Simulate Games";
            this.simulateBtn.UseVisualStyleBackColor = true;
            this.simulateBtn.Click += new System.EventHandler(this.button2_Click);
            // 
            // MainFrm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(580, 435);
            this.Controls.Add(this.simulateBtn);
            this.Controls.Add(this.textBoxResult);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle;
            this.Name = "MainFrm";
            this.Text = "SacBunt";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion
        private System.Windows.Forms.TextBox textBoxResult;
        private System.Windows.Forms.Button simulateBtn;
    }
}

