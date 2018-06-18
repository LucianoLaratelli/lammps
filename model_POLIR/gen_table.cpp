#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

struct program_arguments {
  std::string table_name;
  double theta_count, theta_low, theta_high;
  double r1_count, r1_low, r1_high;
  double r2_count, r2_low, r2_high;
  double theta_const, theta_equil;
  double r_const, r_equil;

};

void print_help()
{
  using std::cout;
  using std::endl;
  cout << "This program takes twelve arguments:" << endl;
  cout << "Table name (string)" << endl;
  cout << "Three angles in degrees: " << endl;
  cout << "Theta count (int) " << endl;
  cout << "Theta low (double) " << endl;
  cout << "Theta high (double) " << endl;
  cout << "Theta equil (double)" << endl;
  cout << "Six distances in angstroms: " << endl;
  cout << "r1 count (int) " << endl;
  cout << "r1 low (double) " << endl;
  cout << "r1 high (double) " << endl;
  cout << "r2 count (int) " << endl;
  cout << "r2 low (double) " << endl;
  cout << "r2 high (double) " << endl;
  cout << "r equil (double) " << endl;
  cout << "Two constants: " << endl;
  cout << "theta const (double) (energy/angle^2)" << endl;
  cout << "r const (double) (energy/angstroms^2)" << endl;

}

program_arguments set_args(std::vector<std::string>args)
{
  if(args.size() != 14) { print_help(); exit(EXIT_FAILURE); }

  program_arguments command_line_args;

  command_line_args.table_name = args.at(0);
  command_line_args.theta_count = std::stoi(args[1]);
  command_line_args.theta_low = std::stod(args[2]) * M_PI / 180.0;
  command_line_args.theta_high = std::stod(args[3]) * M_PI / 180.0;
  command_line_args.theta_equil = std::stod(args[4]) * M_PI / 180.0;
  command_line_args.r1_count = std::stoi(args[5]);
  command_line_args.r1_low = std::stod(args[6]);
  command_line_args.r1_high = std::stod(args[7]);
  command_line_args.r2_count = std::stoi(args[8]);
  command_line_args.r2_low = std::stod(args[9]);
  command_line_args.r2_high = std::stod(args[10]);
  command_line_args.r_equil = std::stod(args[11]);
  command_line_args.theta_const = std::stod(args[12]) * pow((M_PI / 180.0),2);
  command_line_args.r_const = std::stod(args[13]);

  double & theta_c = command_line_args.theta_count;
  double & r1_c = command_line_args.r1_count;
  double & r2_c = command_line_args.r2_count;

  if(theta_c != r1_c && theta_c != r2_c && r1_c != r2_c)
  {
    std::cout << "Sorry at the moment all three variables must have the same count." << std::endl;
    std::cout << "i.e. theta_count = r1_count = r2_count" << std::endl;
    exit(EXIT_FAILURE);
  }

  return command_line_args;
}

double energy(program_arguments user_args, double current_angle, double r1, double r2)
{
  double answer = user_args.theta_const * pow((current_angle-user_args.theta_equil),2);
  answer += user_args.r_const * pow((r1 - user_args.r_equil),2);
  answer += user_args.r_const * pow((r2 - user_args.r_equil),2);
  return answer;
}

double force_theta(program_arguments user_args, double current_quantity)
{
  double answer = 2.0 * user_args.theta_const * (current_quantity - user_args.theta_equil);
  return answer;
}

double force_r(program_arguments user_args, double current_quantity)
{
  double answer = 2.0 * user_args.r_const * (current_quantity - user_args.r_equil);
  return answer;
}
  

int main(int argc, char ** argv)
{
  std::vector<std::string> arguments(argv + 1, argv + argc);

  program_arguments user_args = set_args(arguments);

  std::string table = user_args.table_name;
  transform(table.begin(), table.end(), table.begin(), ::toupper);

  std::ofstream output_file(user_args.table_name+".table");
  output_file << "#angle table" << std::endl << std::endl;
  output_file << table << std::endl;
  output_file << "THETA " << user_args.theta_count << " " << user_args.theta_low << " " << user_args.theta_high;
  output_file << " R1 " << user_args.r1_count << " " << user_args.r1_low << " " << user_args.r1_high;
  output_file << " R2 " << user_args.r2_count << " " << user_args.r2_low << " " << user_args.r2_high;
  output_file << std::endl << std::endl;

  //n angle r1 r2 u f_th f_r1 f_r2

  double dtheta = (user_args.theta_high - user_args.theta_low) / (user_args.theta_count -1);
  double dr1 = (user_args.r1_high - user_args.r1_low) / (user_args.r1_count -1);
  double dr2 = (user_args.r2_high- user_args.r2_low) / (user_args.r2_count -1);

  std::cout << dtheta << " " << dr1 << " " << dr2 << std::endl;

  int count = 0;

  for(double theta = user_args.theta_low; theta <= user_args.theta_high; theta += dtheta)
  {
    for(double r1 = user_args.r1_low; r1 <= user_args.r1_high; r1 += dr1)
    {
      for(double r2 = user_args.r2_low; r2 <= user_args.r2_high; r2 += dr2)
      {
        output_file << count << " " << theta << " "
                    << std::setprecision(6) << r1 << " "
                    << std::setprecision(6) << r2 << " "
                    << std::setprecision(6) << energy(user_args, theta, r1, r2) << " "
                    << std::setprecision(6) << force_theta(user_args, theta) << " "
                    << std::setprecision(6) << force_r(user_args, r1) << " "
                    << std::setprecision(6) << force_r(user_args, r2) << std::endl;
        count++;
      }
    }
  }

  output_file.close();

  int total_count = user_args.theta_count * user_args.r1_count * user_args.r2_count;

  std::cout << "Done! Add the following lines to your LAMMPS input:" << std::endl;
  std::cout << "angle_style table3D linear " <<  total_count << std::endl;
  std::cout << "angle_coeff 1 " << user_args.table_name+".table " << table << std::endl;

  return 0;
}
