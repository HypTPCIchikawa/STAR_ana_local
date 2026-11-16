void draw_canvas_batch(const char* jobsPath);
#ifndef __CINT__
#include "draw_split_canvases.C"

void draw_canvas_batch(const char* jobsPath) {
  try {
    std::ifstream jobs(jobsPath); split_canvas::Need(bool(jobs),"Cannot read PDF job list");
    std::string line; unsigned int count=0;
    while (std::getline(jobs,line)) {
      if (line.empty()) continue;
      const std::vector<std::string> fields=split_canvas::Split(line,'\t');
      split_canvas::Need(fields.size()==5,"Wrong PDF job field count");
      split_canvas::Run(fields[0].c_str(),fields[1].c_str(),fields[2].c_str(),fields[3].c_str(),fields[4].c_str());
      ++count;
    }
    split_canvas::Need(!jobs.bad() && count>0,"Empty/unreadable PDF job list");
    std::cout << "DRAW_CANVAS_BATCH_OK files=" << count << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "DRAW_CANVAS_BATCH_FAILED: " << error.what() << std::endl; gSystem->Exit(2);
  }
}
#endif
