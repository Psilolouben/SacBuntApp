namespace SacBunt
{
    partial class Form1
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
            this.textBoxHomeERA = new System.Windows.Forms.TextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.textBoxHomeRuns = new System.Windows.Forms.TextBox();
            this.textBoxHomeGames = new System.Windows.Forms.TextBox();
            this.textBoxAwayGames = new System.Windows.Forms.TextBox();
            this.textBoxAwayRuns = new System.Windows.Forms.TextBox();
            this.label4 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.textBoxAwayERA = new System.Windows.Forms.TextBox();
            this.textBoxResult = new System.Windows.Forms.TextBox();
            this.button1 = new System.Windows.Forms.Button();
            this.button2 = new System.Windows.Forms.Button();
            this.label7 = new System.Windows.Forms.Label();
            this.textBoxHomeTeam = new System.Windows.Forms.TextBox();
            this.label8 = new System.Windows.Forms.Label();
            this.textBoxAwayTeam = new System.Windows.Forms.TextBox();
            this.SuspendLayout();
            // 
            // textBoxHomeERA
            // 
            this.textBoxHomeERA.Location = new System.Drawing.Point(127, 27);
            this.textBoxHomeERA.Name = "textBoxHomeERA";
            this.textBoxHomeERA.Size = new System.Drawing.Size(100, 20);
            this.textBoxHomeERA.TabIndex = 0;
            this.textBoxHomeERA.Text = "4.18";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(127, 8);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(96, 13);
            this.label1.TabIndex = 1;
            this.label1.Text = "Home Pitcher ERA";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(242, 8);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(93, 13);
            this.label2.TabIndex = 2;
            this.label2.Text = "Home Team Runs";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(353, 8);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(101, 13);
            this.label3.TabIndex = 3;
            this.label3.Text = "Home Team Games";
            // 
            // textBoxHomeRuns
            // 
            this.textBoxHomeRuns.Location = new System.Drawing.Point(245, 27);
            this.textBoxHomeRuns.Name = "textBoxHomeRuns";
            this.textBoxHomeRuns.Size = new System.Drawing.Size(100, 20);
            this.textBoxHomeRuns.TabIndex = 4;
            this.textBoxHomeRuns.Text = "271";
            // 
            // textBoxHomeGames
            // 
            this.textBoxHomeGames.Location = new System.Drawing.Point(356, 27);
            this.textBoxHomeGames.Name = "textBoxHomeGames";
            this.textBoxHomeGames.Size = new System.Drawing.Size(100, 20);
            this.textBoxHomeGames.TabIndex = 5;
            this.textBoxHomeGames.Text = "59";
            // 
            // textBoxAwayGames
            // 
            this.textBoxAwayGames.Location = new System.Drawing.Point(356, 83);
            this.textBoxAwayGames.Name = "textBoxAwayGames";
            this.textBoxAwayGames.Size = new System.Drawing.Size(100, 20);
            this.textBoxAwayGames.TabIndex = 11;
            this.textBoxAwayGames.Text = "64";
            // 
            // textBoxAwayRuns
            // 
            this.textBoxAwayRuns.Location = new System.Drawing.Point(245, 83);
            this.textBoxAwayRuns.Name = "textBoxAwayRuns";
            this.textBoxAwayRuns.Size = new System.Drawing.Size(100, 20);
            this.textBoxAwayRuns.TabIndex = 10;
            this.textBoxAwayRuns.Text = "294";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(353, 64);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(99, 13);
            this.label4.TabIndex = 9;
            this.label4.Text = "Away Team Games";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(242, 64);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(91, 13);
            this.label5.TabIndex = 8;
            this.label5.Text = "Away Team Runs";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(127, 64);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(94, 13);
            this.label6.TabIndex = 7;
            this.label6.Text = "Away Pitcher ERA";
            // 
            // textBoxAwayERA
            // 
            this.textBoxAwayERA.Location = new System.Drawing.Point(127, 83);
            this.textBoxAwayERA.Name = "textBoxAwayERA";
            this.textBoxAwayERA.Size = new System.Drawing.Size(100, 20);
            this.textBoxAwayERA.TabIndex = 6;
            this.textBoxAwayERA.Text = "6.26";
            // 
            // textBoxResult
            // 
            this.textBoxResult.Location = new System.Drawing.Point(13, 180);
            this.textBoxResult.Multiline = true;
            this.textBoxResult.Name = "textBoxResult";
            this.textBoxResult.ScrollBars = System.Windows.Forms.ScrollBars.Vertical;
            this.textBoxResult.Size = new System.Drawing.Size(555, 243);
            this.textBoxResult.TabIndex = 12;
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(15, 122);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(553, 23);
            this.button1.TabIndex = 13;
            this.button1.Text = "Calculate";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(15, 151);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(553, 23);
            this.button2.TabIndex = 14;
            this.button2.Text = "Import From File";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click);
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(10, 8);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(65, 13);
            this.label7.TabIndex = 15;
            this.label7.Text = "Home Team";
            // 
            // textBoxHomeTeam
            // 
            this.textBoxHomeTeam.Location = new System.Drawing.Point(13, 27);
            this.textBoxHomeTeam.Name = "textBoxHomeTeam";
            this.textBoxHomeTeam.Size = new System.Drawing.Size(100, 20);
            this.textBoxHomeTeam.TabIndex = 16;
            this.textBoxHomeTeam.Text = "MIN";
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(12, 64);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(65, 13);
            this.label8.TabIndex = 17;
            this.label8.Text = "Home Team";
            // 
            // textBoxAwayTeam
            // 
            this.textBoxAwayTeam.Location = new System.Drawing.Point(15, 80);
            this.textBoxAwayTeam.Name = "textBoxAwayTeam";
            this.textBoxAwayTeam.Size = new System.Drawing.Size(100, 20);
            this.textBoxAwayTeam.TabIndex = 18;
            this.textBoxAwayTeam.Text = "SEA";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(580, 435);
            this.Controls.Add(this.textBoxAwayTeam);
            this.Controls.Add(this.label8);
            this.Controls.Add(this.textBoxHomeTeam);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.button2);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.textBoxResult);
            this.Controls.Add(this.textBoxAwayGames);
            this.Controls.Add(this.textBoxAwayRuns);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.textBoxAwayERA);
            this.Controls.Add(this.textBoxHomeGames);
            this.Controls.Add(this.textBoxHomeRuns);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.textBoxHomeERA);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle;
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.TextBox textBoxHomeERA;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.TextBox textBoxHomeRuns;
        private System.Windows.Forms.TextBox textBoxHomeGames;
        private System.Windows.Forms.TextBox textBoxAwayGames;
        private System.Windows.Forms.TextBox textBoxAwayRuns;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.TextBox textBoxAwayERA;
        private System.Windows.Forms.TextBox textBoxResult;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Button button2;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.TextBox textBoxHomeTeam;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.TextBox textBoxAwayTeam;
    }
}

