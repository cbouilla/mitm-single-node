#ifndef MITM_UTIL_FILE_SYSTEM
#define MITM_UTIL_FILE_SYSTEM
#include <iostream>
#include <fstream>
#include <filesystem>


/* 0: Folder doens't exist and could NOT be created
 * 1: Folder exists.
 * 2: Folder did not exist but was created during function call.
 */
inline /* since it's defined in a header */
int create_folder_if_not_exist(const std::string& folder_path) {
  if (!std::filesystem::exists(folder_path)) {
    try {
      std::filesystem::create_directories(folder_path);
      std::cout << "Folder created: " << folder_path << std::endl;
      /* we have created a file  */
      return 2;
    } catch (const std::filesystem::filesystem_error& e) {
      std::cerr << "Error creating folder: " << e.what() << std::endl;
      return 0; /* file doesn't exist and we could not create it. */
    }
  } 
  /* File exists: */
  return 1;
}

/* 0: File doens't exist and could NOT be created
 * 1: File exists.
 * 2: File did not exist but was created during function call.
 */
inline /* since it's defined in a header */
int create_file_if_not_exist(const std::string& file_path) {
  if (!std::filesystem::exists(file_path)) {
    std::ofstream outfile(file_path);
    if (outfile.is_open()) {
      std::cout << "File created: " << file_path << std::endl;
      outfile.close();
      return 2;
    } else {
      std::cerr << "Error creating file: Unable to open file." << std::endl;
      return 0;
    }
  }
  return 1; /* File already exist */
}
#endif 
