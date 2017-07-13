using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Net;
using Newtonsoft.Json.Linq;
using System.IO;

namespace SacBunt
{
    public partial class MainFrm : Form
    {
        public MainFrm()
        {
            InitializeComponent();
        }


        private void CalculateValues(double homeERA, double awayERA, int homeRuns, int awayRuns, int homeGames, int awayGames, string homeTeam, string awayTeam, List<int> homeSDEV = null, List<int> awaySDEV = null, List<int> homeHitsStdDev = null, List<int> awayHitsStdDev = null)
        {
            var games = new List<Tuple<double, double>>();
            var roundedGames = new List<Tuple<int, int>>();

            //Number of simulations per game
            var N = 1000000;

            for (int i = 0; i < N; i++)
            {
                //Generate seed
                var seed = System.DateTime.Now.ToFileTime() + i;
                var mySimpleRng = new SimpleRNGNew();
                mySimpleRng.SetSeed((uint)(seed >> 16), (uint)(seed % 4294967296));


                //stddev
                var homeSTDDEV = default(double);
                var awaySTDDEV = default(double);
                var homeHits_STDDEV = default(double);
                var awayHitS_TDDEV = default(double);

                #region " Calculate Standard Deviatons"
                foreach (var h in homeSDEV)
                    homeSTDDEV += Math.Abs(Math.Pow(h - homeERA, 2)) / homeSDEV.Count;

                homeSTDDEV = Math.Sqrt(homeSTDDEV);

                foreach (var h in awaySDEV)
                    awaySTDDEV += Math.Abs(Math.Pow(h - awayERA, 2)) / awaySDEV.Count;

                awaySTDDEV = Math.Sqrt(awaySTDDEV);

                foreach (var h in homeHitsStdDev)
                    homeHits_STDDEV += Math.Abs(Math.Pow(h - ((double)homeRuns / homeGames), 2)) / homeHitsStdDev.Count;

                homeHits_STDDEV = Math.Sqrt(homeHits_STDDEV);

                foreach (var h in awayHitsStdDev)
                    awayHitS_TDDEV += Math.Abs(Math.Pow(h - ((double)awayRuns / awayGames), 2)) / awayHitsStdDev.Count;

                awayHitS_TDDEV = Math.Sqrt(awayHitS_TDDEV);

                #endregion

                #region " Generate Random Numbers"
                double[,] U = new double[4, 1];

                //for (int j = 0; j < 6; j++)
                U[0, 0] = Convert.ToInt32(mySimpleRng.GetNormal((homeRuns / homeGames), homeHits_STDDEV));
                U[1, 0] = Convert.ToInt32(mySimpleRng.GetNormal((awayRuns / awayGames), awayHitS_TDDEV));
                U[2, 0] = Convert.ToInt32(mySimpleRng.GetNormal(homeERA, homeSTDDEV));
                U[3, 0] = Convert.ToInt32(mySimpleRng.GetNormal(awayERA, awaySTDDEV));

                #endregion

                var homeScore = Convert.ToInt32((U[0, 0] + U[3, 0]) / 2);
                var awayScore = Convert.ToInt32((U[1, 0] + U[2, 0]) / 2);

                games.Add(new Tuple<double, double>(Math.Abs(homeScore), Math.Abs(awayScore)));
            }

            #region " Retrieve Stats"
            var isRoundUp = false;
            foreach (var game in games)
            {
                if (isRoundUp)
                    roundedGames.Add(new Tuple<int, int>(Convert.ToInt32(Math.Ceiling(game.Item1)), Convert.ToInt32(Math.Ceiling(game.Item2))));
                else
                    roundedGames.Add(new Tuple<int, int>(Convert.ToInt32(Math.Floor(game.Item1)), Convert.ToInt32(Math.Floor(game.Item2))));
            }

            roundedGames = roundedGames.Where(sd => sd.Item1 != sd.Item2).ToList();

            var homeWins = roundedGames.Count(sd => sd.Item1 > sd.Item2);
            var awayWins = roundedGames.Count(sd => sd.Item1 < sd.Item2);

            var totalRunsPerGame = roundedGames.Select(sd => sd.Item1 + sd.Item2);
            var over7hPerc = ((double)totalRunsPerGame.Count(sd => sd >= 7) / roundedGames.Count) * 100;
            var over75hPerc = ((double)totalRunsPerGame.Count(sd => sd > 7.5) / roundedGames.Count) * 100;
            var over8hPerc = ((double)totalRunsPerGame.Count(sd => sd >= 8) / roundedGames.Count) * 100;
            var over85hPerc = ((double)totalRunsPerGame.Count(sd => sd > 8.5) / roundedGames.Count) * 100;
            var over95hPerc = ((double)totalRunsPerGame.Count(sd => sd > 9.5) / roundedGames.Count) * 100;
            var over9hPerc = ((double)totalRunsPerGame.Count(sd => sd >= 9) / roundedGames.Count) * 100;
            var over105hPerc = ((double)totalRunsPerGame.Count(sd => sd > 10.5) / roundedGames.Count) * 100;
            var over10hPerc = ((double)totalRunsPerGame.Count(sd => sd >= 10) / roundedGames.Count) * 100;
            var over11hPerc = ((double)totalRunsPerGame.Count(sd => sd >= 11) / roundedGames.Count) * 100;
            var over115hPerc = ((double)totalRunsPerGame.Count(sd => sd > 11) / roundedGames.Count) * 100;
            var over12hPerc = ((double)totalRunsPerGame.Count(sd => sd >= 12) / roundedGames.Count) * 100;
            var over125hPerc = ((double)totalRunsPerGame.Count(sd => sd > 12) / roundedGames.Count) * 100;

            var homeWinsPerc = (double)homeWins / roundedGames.Count;
            var awayWinsPerc = (double)awayWins / roundedGames.Count;



            this.textBoxResult.Text += homeTeam + ": " + Math.Round(homeWinsPerc, 2) * 100 + " % " + awayTeam + ": " + Math.Round(awayWinsPerc, 2) * 100 + "% " + Environment.NewLine + "    "
                 + " Over 7: " + Math.Round(over7hPerc, 2) + "% " + Environment.NewLine + "    "
                + " Over 7.5: " + Math.Round(over75hPerc, 2) + "% " + Environment.NewLine + "    "
                + " Over 8: " + Math.Round(over8hPerc, 2) + "% " + Environment.NewLine + "    "
                + " Over 8.5: " + Math.Round(over85hPerc, 2) + "% " + Environment.NewLine + "    "
                + " Over 9: " + Math.Round(over9hPerc, 2) + "% " + Environment.NewLine + "    "
                + "Over 9.5: " + Math.Round(over95hPerc, 2) + "% " + Environment.NewLine + "    "
                + "Over 10: " + Math.Round(over10hPerc, 2) + "% " + Environment.NewLine + "    "
                + "Over 10.5: " + Math.Round(over105hPerc, 2) + "%" + Environment.NewLine + "    "
                + "Over 11: " + Math.Round(over11hPerc, 2) + "%" + Environment.NewLine + "    "
                 + "Over 11.5: " + Math.Round(over115hPerc, 2) + "%" + Environment.NewLine + "    "
                  + "Over 12: " + Math.Round(over125hPerc, 2) + "%" + Environment.NewLine + "    "
                   + "Over 12.5: " + Math.Round(over125hPerc, 2) + "%" + Environment.NewLine + "    ";

            AppendResultToTextFile(homeTeam + ": " + Math.Round(homeWinsPerc, 2) * 100 + " % " + awayTeam + ": " + Math.Round(awayWinsPerc, 2) * 100 + "% " + Environment.NewLine + "    "
                 + " Over 7: " + Math.Round(over7hPerc, 2) + "% " + Environment.NewLine + "    "
                + " Over 7.5: " + Math.Round(over75hPerc, 2) + "% " + Environment.NewLine + "    "
                + " Over 8: " + Math.Round(over8hPerc, 2) + "% " + Environment.NewLine + "    "
                + " Over 8.5: " + Math.Round(over85hPerc, 2) + "% " + Environment.NewLine + "    "
                + " Over 9: " + Math.Round(over9hPerc, 2) + "% " + Environment.NewLine + "    "
                + "Over 9.5: " + Math.Round(over95hPerc, 2) + "% " + Environment.NewLine + "    "
                + "Over 10: " + Math.Round(over10hPerc, 2) + "% " + Environment.NewLine + "    "
                + "Over 10.5: " + Math.Round(over105hPerc, 2) + "%" + Environment.NewLine + "    "
                + "Over 11: " + Math.Round(over11hPerc, 2) + "%" + Environment.NewLine
                 + "Over 11.5: " + Math.Round(over115hPerc, 2) + "%" + Environment.NewLine + "    "
                  + "Over 12: " + Math.Round(over125hPerc, 2) + "%" + Environment.NewLine + "    "
                   + "Over 12.5: " + Math.Round(over125hPerc, 2) + "%" + Environment.NewLine + "    ");

            #endregion
        }


        #region "Helper Methods"
        private List<T> GetStat<T>(Helpers.SearchType type, string id = "", DateTime date = default(DateTime))
        {
            int year, month, day;
            if (date != default(DateTime))
            {
                day = date.Day;
                month = date.Month;
                year = date.Year;

            }
            else
            {
                date = DateTime.Now;
                day = date.Day;
                month = date.Month;
                year = date.Year;
            }

            var daystr = day < 10 ? "0" + day.ToString() : day.ToString();
            var monthstr = month < 10 ? "0" + month.ToString() : month.ToString();

            switch (type)
            {
                case Helpers.SearchType.BaseballGame:
                    var url = @"http://mlb.mlb.com/gdcross/components/game/mlb/year_" + year + "/month_" + monthstr + "/day_" + daystr + "/master_scoreboard.json";
                    var x = new List<Helpers.BaseballGame>();
                    var json = new WebClient().DownloadString(url);

                    JToken token = JObject.Parse(json);

                    var games = token.SelectToken("data.games.game");




                    foreach (var gm in games)
                    {
                        try
                        {
                            var bgame = new Helpers.BaseballGame();
                            bgame.HomeTeam = (string)gm.SelectToken("home_team_city");
                            bgame.AwayTeam = (string)gm.SelectToken("away_team_city");
                            if (gm.SelectToken("home_probable_pitcher.era") != null)
                                bgame.HomePitcherERA = (decimal)gm.SelectToken("home_probable_pitcher.era");
                            if (gm.SelectToken("home_probable_pitcher.id") != null)
                                bgame.HomePitcherID = (int)gm.SelectToken("home_probable_pitcher.id");
                            if (gm.SelectToken("away_probable_pitcher.era") != null)
                                bgame.AwayPitcherERA = (decimal)gm.SelectToken("away_probable_pitcher.era");
                            if (gm.SelectToken("away_probable_pitcher.id") != null)
                                bgame.AwayPitcherID = (int)gm.SelectToken("away_probable_pitcher.id");
                            bgame.HomeTeamID = (int)gm.SelectToken("home_team_id");
                            bgame.AwayTeamID = (int)gm.SelectToken("away_team_id");

                            if (gm.SelectToken("linescore") != null)
                            {
                                bgame.HomeTeamScore = (int)gm.SelectToken("linescore.r.home");
                                bgame.AwayTeamScore = (int)gm.SelectToken("linescore.r.away");
                            }
                            x.Add(bgame);
                        }
                        catch (Exception ex)
                        {
                            continue;
                        }
                    }

                    return (List<T>)Convert.ChangeType(x, typeof(List<T>));
                case Helpers.SearchType.RunsHistory:
                    var x_runs = new List<Helpers.RunsHistory>();
                    url = @"http://mlb.mlb.com/lookup/json/named.sport_pitching_game_log_composed.bam?game_type=%27R%27&league_list_id=%27mlb%27&player_id=" + id + "&season=2017";

                    json = new WebClient().DownloadString(url);

                    token = JObject.Parse(json);

                    var pitchGames = token.SelectToken("sport_pitching_game_log_composed.sport_pitching_game_log.queryResults.row");
                    var p_runs = new Helpers.RunsHistory();
                    p_runs.PitcherID = Convert.ToInt32(id);
                    foreach (var gm in pitchGames)
                    {
                        if (gm.SelectToken("r") != null)
                            p_runs.Runs.Add((int)gm.SelectToken("r"));


                    }
                    x_runs.Add(p_runs);
                    return (List<T>)Convert.ChangeType(x_runs, typeof(List<T>));
                case Helpers.SearchType.TeamRuns:
                    var x_Teamruns = new List<Helpers.TeamRun>();
                    url = @"http://mlb.mlb.com/lookup/json/named.team_hitting_season_leader_master.bam?season=2017&sort_order=%27desc%27&sort_column=%27avg%27&game_type=%27R%27&sport_code=%27mlb%27&recSP=1&recPP=50";
                    json = new WebClient().DownloadString(url);

                    token = JObject.Parse(json);

                    var teamRuns = token.SelectToken("team_hitting_season_leader_master.queryResults.row");

                    foreach (var gm in teamRuns)
                    {
                        var t = new Helpers.TeamRun();

                        if (gm.SelectToken("r") != null)
                        {
                            t.Runs = ((int)gm.SelectToken("r"));
                            t.Games = ((int)gm.SelectToken("g"));
                            t.TeamID = ((int)gm.SelectToken("team_id"));

                            x_Teamruns.Add(t);
                        }


                    }

                    return (List<T>)Convert.ChangeType(x_Teamruns, typeof(List<T>));
                default:
                    return null;

            }


        }

        private void button2_Click(object sender, EventArgs e)
        {
            var path = Application.StartupPath + @"\games.txt";

            using (StreamWriter sw = File.CreateText(path))
            {
                sw.WriteLine(string.Empty);
            }

            #region " Read Web Values "
            var forecast = new Helpers.Forecast();

            forecast.Baseball_Game = this.GetStat<Helpers.BaseballGame>(Helpers.SearchType.BaseballGame, "", DateTime.Now);

            foreach (var phist in forecast.Baseball_Game)
            {
                var pitchingHistoryHome = new List<Helpers.RunsHistory>();
                pitchingHistoryHome = this.GetStat<Helpers.RunsHistory>(Helpers.SearchType.RunsHistory, phist.HomePitcherID.ToString());

                var pitchingHistoryAway = new List<Helpers.RunsHistory>();
                pitchingHistoryAway = this.GetStat<Helpers.RunsHistory>(Helpers.SearchType.RunsHistory, phist.AwayPitcherID.ToString());

                forecast.Runs_History.AddRange(pitchingHistoryHome);
                forecast.Runs_History.AddRange(pitchingHistoryAway);
            }

            forecast.Team_Runs = this.GetStat<Helpers.TeamRun>(Helpers.SearchType.TeamRuns);


            var history = new List<Tuple<int, int>>();
            //History
            var historyCounter = 0;
            var currentDate = DateTime.Now;



            foreach (var game in forecast.Baseball_Game)
            {
                historyCounter = 0;
                currentDate = DateTime.Now;
                while (historyCounter < 10)
                {
                    var datesGames = this.GetStat<Helpers.BaseballGame>(Helpers.SearchType.BaseballGame, "", currentDate.AddDays(-1));

                    var homeResult = datesGames.FirstOrDefault(sd => sd.HomeTeamID == game.HomeTeamID || sd.AwayTeamID == game.HomeTeamID);
                    var awayResult = datesGames.FirstOrDefault(sd => sd.HomeTeamID == game.AwayTeamID || sd.AwayTeamID == game.AwayTeamID);
                    if (homeResult != null)
                    {
                        var homeScore = game.HomeTeamID == homeResult.HomeTeamID ? homeResult.HomeTeamScore : homeResult.AwayTeamScore;
                        history.Add(new Tuple<int, int>(game.HomeTeamID, homeScore));
                    }
                    if (awayResult != null)
                    {
                        var awayScore = game.AwayTeamID == awayResult.HomeTeamID ? awayResult.HomeTeamScore : awayResult.AwayTeamScore;
                        history.Add(new Tuple<int, int>(game.AwayTeamID, awayScore));
                    }

                    historyCounter++;
                    currentDate = currentDate.AddDays(-1);
                }
            }

            #endregion

            foreach (var game in forecast.Baseball_Game)
            {
                var homeRuns = forecast.Team_Runs.FirstOrDefault(sd => sd.TeamID == game.HomeTeamID).Runs;
                var awayRuns = forecast.Team_Runs.FirstOrDefault(sd => sd.TeamID == game.AwayTeamID).Runs;
                var homeGames = forecast.Team_Runs.FirstOrDefault(sd => sd.TeamID == game.HomeTeamID).Games;
                var awayGames = forecast.Team_Runs.FirstOrDefault(sd => sd.TeamID == game.AwayTeamID).Games;

                var homePitcherRuns = forecast.Runs_History.FirstOrDefault(sd => sd.PitcherID == game.HomePitcherID).Runs.ToList();
                var awayPitcherRuns = forecast.Runs_History.FirstOrDefault(sd => sd.PitcherID == game.AwayPitcherID).Runs.ToList();

                var homeTeamRunsHistory = history.Where(sd => sd.Item1 == game.HomeTeamID).Select(sd => sd.Item2);
                var awayTeamRunsHistory = history.Where(sd => sd.Item1 == game.AwayTeamID).Select(sd => sd.Item2);

                if (homePitcherRuns.Count() > 0 && awayPitcherRuns.Count() > 0)
                    CalculateValues(Convert.ToDouble(game.HomePitcherERA), Convert.ToDouble(game.AwayPitcherERA), homeRuns, awayRuns, homeGames, awayGames, game.HomeTeam, game.AwayTeam, homePitcherRuns, awayPitcherRuns, homeTeamRunsHistory.ToList(), awayTeamRunsHistory.ToList());
            }

            //var path = Application.StartupPath + @"\games.txt";

            //string line;

            //System.IO.StreamReader file = new System.IO.StreamReader(path);
            //while ((line = file.ReadLine()) != null)
            //{
            //    var args = line.Split('/');
            //    var homeTeam = args[0];
            //    var awayTeam = args[1];
            //    var homeERA = Convert.ToDouble(args[2]);
            //    var homeRuns = Convert.ToInt32(args[3]);
            //    var homeGames = Convert.ToInt32(args[4]);
            //    var awayERA = Convert.ToDouble(args[5]);
            //    var awayRuns = Convert.ToInt32(args[6]);
            //    var awayGames = Convert.ToInt32(args[7]);
            //    var homeStdDev = args[8].Split('-').Select(sf => Convert.ToInt32(sf)).ToList();
            //    var awayStdDev = args[9].Split('-').Select(sf => Convert.ToInt32(sf)).ToList();
            //    var homeHitsStdDev = args[10].Split('-').Select(sf => Convert.ToInt32(sf)).ToList();
            //    var awayHitsStdDev = args[11].Split('-').Select(sf => Convert.ToInt32(sf)).ToList();

            //    CalculateValues(homeERA, awayERA, homeRuns, awayRuns, homeGames, awayGames, homeTeam, awayTeam, homeStdDev, awayStdDev, homeHitsStdDev, awayHitsStdDev);

            //}

            //file.Close();

        }

        private void AppendResultToTextFile(string result)
        {

            var path = Application.StartupPath + @"\games.txt";


            using (StreamWriter sw = File.AppendText(path))
            {
                sw.WriteLine(result);
            }



        }
        #endregion

    }
}