using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SacBunt
{
    public class Helpers
    {
        public enum SearchType
        {
            BaseballGame = 1,
            RunsHistory = 2,
            TeamRuns = 3
        }

        public class Forecast
        {
            public Forecast()
            {
                Baseball_Game = new List<BaseballGame>();
                Runs_History = new List<RunsHistory>();
                Team_Runs = new List<TeamRun>();
            }

            public List<BaseballGame> Baseball_Game { get; set; }
            public List<RunsHistory> Runs_History { get; set; }
            public List<TeamRun> Team_Runs { get; set; }

        }

        public class BaseballGame
        {
            public int ID { get; set; }
            public string HomeTeam { get; set; }
            public string AwayTeam { get; set; }
            public decimal HomePitcherERA { get; set; }
            public decimal AwayPitcherERA { get; set; }
            public int HomePitcherID { get; set; }
            public int AwayPitcherID { get; set; }
            public int HomeTeamID { get; set; }
            public int AwayTeamID { get; set; }
            public int HomeTeamScore { get; set; }
            public int AwayTeamScore { get; set; }
        }

        public class RunsHistory
        {
            public RunsHistory()
            {
                Runs = new List<int>();
            }

            public int PitcherID { get; set; }
            public List<int> Runs { get; set; }

        }

        public class TeamRun
        {
            public int TeamID { get; set; }
            public int Runs { get; set; }
            public int Games { get; set;}
        }
    }
}
